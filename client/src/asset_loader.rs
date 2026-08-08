use std::{
    path::{Path, PathBuf},
    rc::Rc,
    sync::Arc,
    thread::{self, JoinHandle},
};

use ash::vk;
use lahar::{ParallelQueue, parallel_queue};
use skid_steer::Context;
use tokio_util::sync::CancellationToken;

use crate::{
    Config,
    graphics::Base,
    growable_ring::{self, GrowableRing},
};

#[derive(Clone, Copy)]
struct HelperSemaphore(vk::Semaphore);

impl HelperSemaphore {
    unsafe fn new(device: &ash::Device) -> Self {
        HelperSemaphore(unsafe {
            device
                .create_semaphore(
                    &vk::SemaphoreCreateInfo::default().push_next(
                        &mut vk::SemaphoreTypeCreateInfo::default()
                            .semaphore_type(vk::SemaphoreType::TIMELINE)
                            .initial_value(0),
                    ),
                    None,
                )
                .unwrap()
        })
    }

    unsafe fn get_next_value(&self, device: &ash::Device) -> u64 {
        unsafe { device.get_semaphore_counter_value(self.0).unwrap() + 1 }
    }

    unsafe fn signal(&self, device: &ash::Device) {
        let next_value = unsafe { self.get_next_value(device) };
        unsafe {
            device
                .signal_semaphore(
                    &vk::SemaphoreSignalInfo::default()
                        .semaphore(self.0)
                        .value(next_value),
                )
                .unwrap()
        }
    }

    unsafe fn destroy(self, device: &ash::Device) {
        unsafe { device.destroy_semaphore(self.0, None) };
    }
}

/// Contains all the dependencies necessary to load assets.
pub struct AssetLoadContext {
    gfx: Arc<Base>,
    config: Arc<Config>,
    parallel_queue_handle: parallel_queue::Handle,
    parallel_queue_semaphore_watcher: tokio::sync::watch::Receiver<u64>,
    queue_unpark_semaphore: HelperSemaphore,
    staging: Arc<GrowableRing>,
}

// Rather than exposing its fields, we expose helper functions for the kinds of tasks
// one would need the context for. This helps add some level of separation between how global
// data is organized and what asset loading code sees
impl AssetLoadContext {
    /// # Safety
    /// - [`Work::cmd`] must not be used outside the lifetime of the returned [`Work`]
    /// - Any Vulkan resources this work uses must not be destroyed before the [`Work`] is dropped
    pub unsafe fn begin_work(&self) -> parallel_queue::Work<'_> {
        unsafe { self.parallel_queue_handle.begin(&self.gfx.device) }
    }

    pub fn alloc<T>(
        &self,
        count: usize,
        align: usize,
        free_at: u64,
    ) -> growable_ring::Allocation<T> {
        self.staging.alloc(&self.gfx, count, align, free_at)
    }

    pub fn device(&self) -> &ash::Device {
        self.gfx.device.as_ref()
    }

    pub fn memory_properties(&self) -> &vk::PhysicalDeviceMemoryProperties {
        &self.gfx.memory_properties
    }

    pub fn queue_family(&self) -> u32 {
        self.gfx.queue_family
    }

    pub async fn wait_for_completion(&self, semaphore_value: u64) {
        // To actually get the work to start, we need to unpark the queue. We do it here to avoid getting stuck awaiting something we never kicked off.
        unsafe { self.queue_unpark_semaphore.signal(&self.gfx.device) };
        if self
            .parallel_queue_semaphore_watcher
            .clone()
            .wait_for(|&value| value >= semaphore_value)
            .await
            .is_err()
        {
            // Don't assume things have finished loading just because `parallel_queue_semaphore_notifier` was dropped.
            // Instead, print an error and wait forever, since the work will never complete.
            tracing::error!(
                "Work was submitted to the parallel queue but never completed. This should never happen."
            );
            std::future::pending::<()>().await;
        }
    }

    pub fn block_on_work_completion(&self, semaphore_value: u64) {
        // To actually get the work to start, we need to unpark the queue. We do it here to avoid getting stuck awaiting something we never kicked off.
        unsafe { self.queue_unpark_semaphore.signal(&self.gfx.device) };

        // Blocking on an async function is tricky, as doing it naively can cause Tokio to panic if we run into a situation where
        // a runtime is used within another runtime. Because of that, rather than using `parallel_queue_semaphore_watcher`,
        // we just use the timeline semaphore directly.
        unsafe {
            self.gfx
                .device
                .wait_semaphores(
                    &vk::SemaphoreWaitInfo::default()
                        .semaphores(&[self.parallel_queue_handle.semaphore()])
                        .values(&[semaphore_value]),
                    !0,
                )
                .unwrap()
        };
    }

    pub fn find_asset(&self, path: &Path) -> Option<PathBuf> {
        self.config.find_asset(path)
    }
}

pub struct AssetLoader {
    gfx: Arc<Base>,
    shutdown_token: CancellationToken,
    queue_unpark_semaphore: HelperSemaphore,
    loader: skid_steer::Loader,
    task_executor_threads: Vec<JoinHandle<parallel_queue::Handle>>,
    parallel_queue_driver_thread: Option<JoinHandle<ParallelQueue>>,
    parallel_queue_semaphore_notifier: tokio::sync::watch::Sender<u64>,
    staging: Arc<GrowableRing>,
}

impl AssetLoader {
    pub fn new(gfx: Arc<Base>, config: Arc<Config>) -> Self {
        let shutdown_token = CancellationToken::new();
        let queue = unsafe { ParallelQueue::new(&gfx.device, gfx.queue_family, gfx.queue, None) };
        let loader = skid_steer::Loader::new();
        let staging = Arc::new(GrowableRing::new(&gfx, 32 * 1024 * 1024));
        let queue_unpark_semaphore = unsafe { HelperSemaphore::new(&gfx.device) };
        let (parallel_queue_semaphore_notifier, parallel_queue_semaphore_watcher) =
            tokio::sync::watch::channel(0);

        let mut task_executor_threads = vec![];
        tracing::debug!(
            "Using asset load parallelism {}",
            config.asset_load_parallelism
        );
        for i in 0..(config.asset_load_parallelism) {
            let asset_load_context = AssetLoadContext {
                gfx: Arc::clone(&gfx),
                config: Arc::clone(&config),
                parallel_queue_handle: unsafe { queue.handle(&gfx.device) },
                parallel_queue_semaphore_watcher: parallel_queue_semaphore_watcher.clone(),
                queue_unpark_semaphore,
                staging: Arc::clone(&staging),
            };
            let loader = loader.clone();

            let thread = thread::Builder::new()
                .name(format!("task_executor_{}", i).to_owned())
                .spawn(move || run_task_executor(asset_load_context, loader))
                .unwrap();
            task_executor_threads.push(thread);
        }

        let parallel_queue_driver_thread = Some({
            let gfx = Arc::clone(&gfx);
            let staging = Arc::clone(&staging);
            let shutdown_token = shutdown_token.clone();
            let parallel_queue_semaphore_notifier = parallel_queue_semaphore_notifier.clone();

            thread::Builder::new()
                .name("parallel_queue_driver".to_owned())
                .spawn(move || {
                    run_parallel_queue_driver(
                        &gfx,
                        queue,
                        queue_unpark_semaphore,
                        shutdown_token,
                        &parallel_queue_semaphore_notifier,
                        &staging,
                    )
                })
                .unwrap()
        });

        AssetLoader {
            gfx,
            shutdown_token,
            queue_unpark_semaphore,
            loader,
            task_executor_threads,
            parallel_queue_driver_thread,
            parallel_queue_semaphore_notifier,
            staging,
        }
    }

    pub fn load<S: skid_steer::Source>(
        &self,
        source: S,
    ) -> skid_steer::Asset<<S as skid_steer::Source>::Output> {
        self.loader.load(source)
    }
}

#[must_use]
fn run_task_executor(
    asset_load_context: AssetLoadContext,
    loader: skid_steer::Loader,
) -> parallel_queue::Handle {
    let asset_load_context = Rc::new(asset_load_context);

    tokio::runtime::LocalRuntime::new()
        .unwrap()
        .block_on(async {
            while let Some(task) = loader.next_task().await {
                tracing::trace!(
                    "Found task on {}",
                    thread::current().name().unwrap_or("<unnamed>")
                );
                let asset_load_context = Rc::clone(&asset_load_context);
                tokio::task::spawn_local(async move {
                    let mut context = Context::new();
                    context.insert::<AssetLoadContext>(&asset_load_context);
                    task.run(&context).await;
                    tracing::trace!(
                        "Task complete on {}",
                        thread::current().name().unwrap_or("<unnamed>")
                    );
                });
            }
            tracing::trace!(
                "Ending task executor {}",
                thread::current().name().unwrap_or("<unnamed>")
            );
        });
    let asset_load_context = Rc::try_unwrap(asset_load_context)
        .map_err(|_| ())
        .expect("runtime using this context should already be dropped");
    asset_load_context.parallel_queue_handle
}

fn run_parallel_queue_driver(
    gfx: &Base,
    mut queue: ParallelQueue,
    queue_unpark_semaphore: HelperSemaphore,
    shutdown_token: CancellationToken,
    parallel_queue_semaphore_notifier: &tokio::sync::watch::Sender<u64>,
    staging: &GrowableRing,
) -> ParallelQueue {
    loop {
        // Systems increment the `queue_unpark_semaphore` value when they want to guarantee
        // that we don't park unless certain things are done. `ParallelQueueWaiter`
        // wants to ensure that `ParallelQueue::drive` is called, while the cleanup code
        // wants to ensure `cancellation_recieve` is checked. Therefore, we put these two
        // operations between reading and waiting on  `queue_unpark_semaphore`
        let queue_unpark_semaphore_next_value =
            unsafe { queue_unpark_semaphore.get_next_value(&gfx.device) };
        unsafe { queue.drive(&gfx.device) };
        if shutdown_token.is_cancelled() {
            break;
        }
        let semaphore_value = unsafe {
            queue.park(
                &gfx.device,
                queue_unpark_semaphore.0,
                queue_unpark_semaphore_next_value,
            )
        };
        let _ = parallel_queue_semaphore_notifier.send(semaphore_value);
        unsafe { staging.tick(&gfx.device, semaphore_value) };
    }
    queue
}

impl Drop for AssetLoader {
    fn drop(&mut self) {
        self.shutdown_token.cancel();
        unsafe { self.queue_unpark_semaphore.signal(&self.gfx.device) };
        let mut queue = self
            .parallel_queue_driver_thread
            .take()
            .unwrap()
            .join()
            .unwrap();
        tracing::trace!("Shutting down down AssetLoader");
        self.loader.close();
        let parallel_queue_handles: Vec<_> = self
            .task_executor_threads
            .drain(..)
            .map(|thread| thread.join().unwrap())
            .collect();
        unsafe { queue.drain(&self.gfx.device) };
        // We need to notify for the final value the parallel queue semaphore reaches.
        let _ = self.parallel_queue_semaphore_notifier.send(unsafe {
            self.gfx
                .device
                .get_semaphore_counter_value(queue.semaphore())
                .unwrap()
        });
        for mut handle in parallel_queue_handles {
            // Safety: The queue has been drained, and no handles are left behind, so no work should be able to be in flight.
            unsafe { handle.destroy(&self.gfx.device) };
        }
        unsafe { queue.destroy(&self.gfx.device) };
        unsafe { self.queue_unpark_semaphore.destroy(&self.gfx.device) };
        unsafe {
            Arc::get_mut(&mut self.staging)
                .expect("All threads using staging should now be joined")
                .destroy(&self.gfx.device)
        };
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashSet, sync::Mutex, time::Duration};

    use super::*;

    /// Record-keeping for what has happened so far with the asset loader, useful for assertions
    #[derive(Debug, Clone)]
    struct EventList {
        events: Arc<Mutex<Vec<Event>>>,
        num_checked_events: usize,
    }

    impl EventList {
        fn new() -> Self {
            EventList {
                events: Arc::new(Mutex::new(vec![])),
                num_checked_events: 0,
            }
        }

        fn push_event(&self, event: Event) {
            self.events.lock().unwrap().push(event);
        }

        fn get_all_events(&mut self) -> Vec<Event> {
            let events = self.events.lock().unwrap().clone();
            self.num_checked_events = events.len();
            events
        }

        /// Drains all events that were returned by the most recent call to get_all_events
        fn drain_queried_events(&mut self) {
            self.events
                .lock()
                .unwrap()
                .drain(0..self.num_checked_events);
            self.num_checked_events = 0;
        }
    }

    #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
    enum Event {
        Progress {
            asset: String,
            percent_progress: u32,
        },
        Loaded {
            asset: String,
        },
        Freed {
            asset: String,
        },
        LoadCanceled {
            asset: String,
        },
    }

    impl Event {
        fn progress(asset: &str, percent_progress: u32) -> Self {
            Event::Progress {
                asset: asset.to_owned(),
                percent_progress,
            }
        }

        fn loaded(asset: &str) -> Self {
            Event::Loaded {
                asset: asset.to_owned(),
            }
        }

        fn freed(asset: &str) -> Self {
            Event::Freed {
                asset: asset.to_owned(),
            }
        }

        fn load_cancelled(asset: &str) -> Self {
            Event::LoadCanceled {
                asset: asset.to_owned(),
            }
        }
    }

    struct TestAsset<T: Send + 'static> {
        asset: skid_steer::Asset<T>,
        progress_sender: tokio::sync::mpsc::UnboundedSender<u32>,
    }

    impl<T: Send + 'static> TestAsset<T> {
        fn add_percent_progress(&self, percent_progress: u32) {
            self.progress_sender.send(percent_progress).unwrap();
        }

        fn wait_for_completion(&self) {
            // Note: This design pattern is somewhat dangerous because it can panic if called from within another runtime.
            // Since this is only used in unit tests, this should be an acceptable level of risk.
            tokio::runtime::LocalRuntime::new()
                .unwrap()
                .block_on(async {
                    tokio::select! {
                        biased;
                        _ = self.asset.get() => (),
                        _ = tokio::time::sleep(Duration::from_secs(5)) => { panic!("Timed out waiting for asset to load") },
                    };
                });
        }
    }

    struct DummyAsset {
        name: String,
        events: EventList,
    }

    struct DummyAssetSource {
        name: String,
        events: EventList,
        progress_receiver: tokio::sync::mpsc::UnboundedReceiver<u32>,
    }

    impl skid_steer::Source for DummyAssetSource {
        type Output = DummyAsset;

        async fn load(mut self, context: &Context<'_>) -> Option<Self::Output> {
            let mut status = DummyAssetStatus {
                name: self.name.clone(),
                progress: 0,
                can_cancel: true,
                events: self.events.clone(),
            };

            // Use the parallel queue and set up an allocation to exercise this functionality
            let ctx: &AssetLoadContext = context.get().unwrap();
            let work = unsafe { ctx.begin_work() };
            let finish_time = work.time().get();
            let _alloc = ctx.alloc::<u8>(8, 1, finish_time);

            while status.progress < 100 {
                status.progress += self.progress_receiver.recv().await?;
                self.events
                    .push_event(Event::progress(&self.name, status.progress));
            }

            work.end();
            ctx.wait_for_completion(finish_time).await;

            status.can_cancel = false; // Done loading
            self.events.push_event(Event::loaded(&self.name));
            Some(DummyAsset {
                name: self.name,
                events: self.events,
            })
        }

        fn free(output: Self::Output, _context: &Context) {
            output.events.push_event(Event::freed(&output.name));
        }
    }

    struct DummyAssetStatus {
        name: String,
        progress: u32,
        can_cancel: bool,
        events: EventList,
    }

    impl Drop for DummyAssetStatus {
        fn drop(&mut self) {
            if self.can_cancel {
                self.events.push_event(Event::load_cancelled(&self.name));
            }
        }
    }

    fn init_asset_loader(asset_load_parallelism: u32) -> AssetLoader {
        let gfx = Arc::new(Base::headless());
        let config = Arc::new({
            let mut config = Config::create_for_test();
            config.asset_load_parallelism = asset_load_parallelism;
            config
        });
        AssetLoader::new(Arc::clone(&gfx), Arc::clone(&config))
    }

    fn load_dummy_asset(
        asset_loader: &AssetLoader,
        events: &EventList,
        name: &str,
    ) -> TestAsset<DummyAsset> {
        let (progress_sender, progress_receiver) = tokio::sync::mpsc::unbounded_channel();
        TestAsset {
            asset: asset_loader.load(DummyAssetSource {
                name: name.to_owned(),
                events: events.clone(),
                progress_receiver,
            }),
            progress_sender,
        }
    }

    fn test_eventual_success(mut f: impl FnMut() -> Result<(), anyhow::Error>) {
        let timeout = Duration::from_secs(3);
        if f().is_ok() {
            return;
        }
        let start_time = std::time::Instant::now();
        let mut next_poll = Duration::from_millis(16);
        while next_poll < timeout {
            if let Some(sleep_time) = next_poll.checked_sub(start_time.elapsed()) {
                thread::sleep(sleep_time);
            }
            if f().is_ok() {
                return;
            }
            next_poll *= 2;
        }
        if let Some(sleep_time) = timeout.checked_sub(start_time.elapsed()) {
            thread::sleep(sleep_time);
        }
        if let Err(e) = f() {
            panic!("{}", e);
        }
    }

    #[test]
    fn test_load_and_free() {
        let mut events = EventList::new();
        let asset_loader = init_asset_loader(2);
        let dummy_asset = load_dummy_asset(&asset_loader, &events, "asset");
        dummy_asset.add_percent_progress(50);
        dummy_asset.add_percent_progress(50);
        dummy_asset.wait_for_completion();
        assert!(dummy_asset.asset.try_get().is_some());
        assert_eq!(
            events.get_all_events(),
            &[
                Event::progress("asset", 50),
                Event::progress("asset", 100),
                Event::loaded("asset")
            ]
        );
        events.drain_queried_events();
        drop(dummy_asset);
        drop(asset_loader);
        assert_eq!(events.get_all_events(), &[Event::freed("asset")]);
    }

    #[test]
    fn test_concurrency_and_cancellation() {
        let mut events = EventList::new();
        let asset_loader = init_asset_loader(2);
        let assets: Vec<_> = (0..4)
            .map(|i| load_dummy_asset(&asset_loader, &events, &format!("asset{i}")))
            .collect();
        for asset in &assets {
            asset.add_percent_progress(50);
        }
        test_eventual_success(|| {
            let actual_events = events.get_all_events().into_iter().collect::<HashSet<_>>();
            let expected_events = [
                Event::progress("asset0", 50),
                Event::progress("asset1", 50),
                Event::progress("asset2", 50),
                Event::progress("asset3", 50),
            ]
            .into_iter()
            .collect::<HashSet<_>>();
            anyhow::ensure!(
                actual_events == expected_events,
                "Expected: {expected_events:?}\nActual: {actual_events:?}"
            );
            Ok(())
        });
    }
}
