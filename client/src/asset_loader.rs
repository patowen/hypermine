use std::{
    path::{Path, PathBuf},
    rc::Rc,
    sync::Arc,
    thread::{self, JoinHandle},
};

use ash::vk;
use lahar::{ParallelQueue, parallel_queue};
use skid_steer::Context;

use crate::{
    Config,
    graphics::Base,
    growable_ring::{self, GrowableRing},
};

#[derive(Clone)]
pub struct ParallelQueueWaiter {
    semaphore: vk::Semaphore,
    queue_unpark_semaphore: HelperSemaphore,
    progress_changed: Arc<tokio::sync::Notify>,
}

impl ParallelQueueWaiter {
    pub async fn wait_for_semaphore(&self, device: &ash::Device, value: u64) {
        // Since wait_for_semaphore is called after we submitted something to the queue,
        // we want the `drive` function of `ParallelQueue` to see the new data, so we should
        // unpark it.
        unsafe { self.queue_unpark_semaphore.signal(device) };
        loop {
            let notified = self.progress_changed.notified();
            let current_progress =
                unsafe { device.get_semaphore_counter_value(self.semaphore).unwrap() };
            if current_progress >= value {
                return;
            }
            notified.await;
        }
    }
}

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
    parallel_queue_waiter: ParallelQueueWaiter,
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
        free_at_completion: &parallel_queue::Work<'_>,
    ) -> growable_ring::Allocation<T> {
        self.staging
            .alloc(&self.gfx, count, align, free_at_completion.time().get())
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

    /// Submits the given work immediately and returns a future you can await to ensure the work completes
    pub fn complete_work(&self, work: parallel_queue::Work<'_>) -> impl Future {
        let finish_time = work.time().get();
        work.end();
        self.parallel_queue_waiter
            .wait_for_semaphore(&self.gfx.device, finish_time)
    }

    pub fn find_asset(&self, path: &Path) -> Option<PathBuf> {
        self.config.find_asset(path)
    }
}

pub struct AssetLoader {
    gfx: Arc<Base>,
    cancellation_send: std::sync::mpsc::Sender<()>, // TODO: Use CancellationToken
    queue_unpark_semaphore: HelperSemaphore,
    loader: skid_steer::Loader,
    task_executor_threads: Vec<JoinHandle<()>>,
    parallel_queue_driver_thread: Option<JoinHandle<ParallelQueue>>,
    staging: Arc<GrowableRing>,
}

impl AssetLoader {
    pub fn new(gfx: Arc<Base>, config: Arc<Config>) -> Self {
        let (cancellation_send, cancellation_receive) = std::sync::mpsc::channel::<()>();
        let queue = unsafe { ParallelQueue::new(&gfx.device, gfx.queue_family, gfx.queue, None) };
        let loader = skid_steer::Loader::new();
        let staging = Arc::new(GrowableRing::new(&gfx, 32 * 1024 * 1024));
        let queue_unpark_semaphore = unsafe { HelperSemaphore::new(&gfx.device) };

        let parallel_queue_waiter = ParallelQueueWaiter {
            semaphore: queue.semaphore(),
            queue_unpark_semaphore,
            progress_changed: Arc::new(tokio::sync::Notify::new()),
        };

        // TODO: Use dynamic number of threads
        let mut task_executor_threads = vec![];
        tracing::debug!(
            "Using asset load parallelism {}",
            config.asset_load_parallelism
        );
        for i in 0..(config.asset_load_parallelism) {
            let gfx = Arc::clone(&gfx);
            let config = Arc::clone(&config);
            let handle = unsafe { queue.handle(&gfx.device) };
            let staging = Arc::clone(&staging);
            let parallel_queue_waiter = parallel_queue_waiter.clone();
            let loader = loader.clone();

            let thread = thread::Builder::new()
                .name(format!("task_executor_{}", i).to_owned())
                .spawn(move || {
                    run_task_executor(gfx, config, handle, staging, parallel_queue_waiter, loader);
                })
                .unwrap();
            task_executor_threads.push(thread);
        }

        let parallel_queue_driver_thread = Some({
            let gfx = Arc::clone(&gfx);
            let staging = Arc::clone(&staging);

            thread::Builder::new()
                .name("parallel_queue_driver".to_owned())
                .spawn(move || {
                    run_parallel_queue_driver(
                        &gfx,
                        queue,
                        queue_unpark_semaphore,
                        cancellation_receive,
                        &parallel_queue_waiter,
                        &staging,
                    )
                })
                .unwrap()
        });

        AssetLoader {
            gfx,
            cancellation_send,
            queue_unpark_semaphore,
            loader,
            task_executor_threads,
            parallel_queue_driver_thread,
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

fn run_task_executor(
    gfx: Arc<Base>,
    config: Arc<Config>,
    handle: parallel_queue::Handle,
    staging: Arc<GrowableRing>,
    parallel_queue_waiter: ParallelQueueWaiter,
    loader: skid_steer::Loader,
) {
    let mut asset_load_context = Rc::new(AssetLoadContext {
        gfx,
        config,
        parallel_queue_handle: handle,
        parallel_queue_waiter,
        staging,
    });
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
    let asset_load_context = Rc::get_mut(&mut asset_load_context)
        .expect("runtime using this context should already be dropped");
    unsafe {
        asset_load_context
            .parallel_queue_handle
            .destroy(&asset_load_context.gfx.device)
    };
}

fn run_parallel_queue_driver(
    gfx: &Base,
    mut queue: ParallelQueue,
    queue_unpark_semaphore: HelperSemaphore,
    cancellation_receive: std::sync::mpsc::Receiver<()>,
    parallel_queue_waiter: &ParallelQueueWaiter,
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
        if cancellation_receive.try_recv() != Err(std::sync::mpsc::TryRecvError::Empty) {
            break;
        }
        let semaphore_value = unsafe {
            queue.park(
                &gfx.device,
                queue_unpark_semaphore.0,
                queue_unpark_semaphore_next_value,
            )
        };
        parallel_queue_waiter.progress_changed.notify_waiters();
        unsafe { staging.tick(&gfx.device, semaphore_value) };
    }
    queue
}

impl Drop for AssetLoader {
    fn drop(&mut self) {
        let _ = self.cancellation_send.send(());
        unsafe { self.queue_unpark_semaphore.signal(&self.gfx.device) };
        let mut queue = self
            .parallel_queue_driver_thread
            .take()
            .unwrap()
            .join()
            .unwrap();
        tracing::trace!("Shutting down down AssetLoader");
        self.loader.close();
        for task_executor_thread in self.task_executor_threads.drain(..) {
            task_executor_thread.join().unwrap();
        }
        unsafe { queue.drain(&self.gfx.device) };
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
    struct EventList(Arc<Mutex<Vec<Event>>>);

    impl EventList {
        fn new() -> Self {
            EventList(Arc::new(Mutex::new(vec![])))
        }

        fn push_event(&self, event: Event) {
            self.0.lock().unwrap().push(event);
        }

        fn drain(&self) -> Vec<Event> {
            self.0.lock().unwrap().drain(..).collect()
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

        async fn load(mut self, _context: &Context<'_>) -> Option<Self::Output> {
            let mut status = DummyAssetStatus {
                name: self.name.clone(),
                progress: 0,
                can_cancel: true,
                events: self.events.clone(),
            };
            while status.progress < 100 {
                status.progress += self.progress_receiver.recv().await?;
                self.events
                    .push_event(Event::progress(&self.name, status.progress));
            }
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

    #[test]
    fn sample_test() {
        let events = EventList::new();
        let asset_loader = init_asset_loader(2);
        let dummy_asset = load_dummy_asset(&asset_loader, &events, "asset");
        dummy_asset.add_percent_progress(50);
        dummy_asset.add_percent_progress(50);
        dummy_asset.wait_for_completion();
        assert!(dummy_asset.asset.try_get().is_some());
        assert_eq!(
            events.drain(),
            &[
                Event::progress("asset", 50),
                Event::progress("asset", 100),
                Event::loaded("asset")
            ]
        );
        drop(dummy_asset);
        drop(asset_loader);
        assert_eq!(events.drain(), &[Event::freed("asset")]);
    }

    #[test]
    fn test_concurrency() {
        let events = EventList::new();
        let asset_loader = init_asset_loader(2);
        let assets: Vec<_> = (0..4)
            .map(|i| load_dummy_asset(&asset_loader, &events, &format!("asset{i}")))
            .collect();
        for asset in &assets {
            asset.add_percent_progress(50);
        }
        std::thread::sleep(Duration::from_secs(1)); // TODO: Use a helper function for the 1-second timeout
        assert_eq!(
            events.drain().into_iter().collect::<HashSet<_>>(),
            [
                Event::progress("asset0", 50),
                Event::progress("asset1", 50),
                Event::progress("asset2", 50),
                Event::progress("asset3", 50)
            ]
            .into_iter()
            .collect::<HashSet<_>>()
        );
    }
}
