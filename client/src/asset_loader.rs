use std::{
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

pub struct AssetLoadContext {
    gfx: Arc<Base>,
    parallel_queue_handle: parallel_queue::Handle,
    parallel_queue_waiter: ParallelQueueWaiter,
    staging: Arc<GrowableRing>,
}

impl AssetLoadContext {
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

    /// Submits the given work immediately and returns a future you can await to ensure the work completes
    pub fn complete_work(&self, work: parallel_queue::Work<'_>) -> impl Future {
        let finish_time = work.time().get();
        work.end();
        self.parallel_queue_waiter
            .wait_for_semaphore(&self.gfx.device, finish_time)
    }
}

pub struct AssetLoader {
    gfx: Arc<Base>,
    cancellation_send: std::sync::mpsc::Sender<()>, // TODO: Use CancellationToken
    queue_unpark_semaphore: HelperSemaphore,
    loader: skid_steer::Loader,
    task_executor_threads: Vec<JoinHandle<()>>,
    parallel_queue_driver_thread: Option<JoinHandle<ParallelQueue>>,
    staging: Option<Arc<GrowableRing>>,
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
        for i in 0..2 {
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
            staging: Some(staging),
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
    let runtime = tokio::runtime::LocalRuntime::new().unwrap();
    let asset_load_context = Rc::new(AssetLoadContext {
        gfx,
        parallel_queue_handle: handle,
        parallel_queue_waiter,
        staging,
    });
    runtime.block_on(async {
        while let Some(task) = loader.next_task().await {
            tracing::trace!(
                "Found task on {}",
                thread::current().name().unwrap_or("<unnamed>")
            );
            let config = Arc::clone(&config);
            let asset_load_context = Rc::clone(&asset_load_context);
            tokio::task::spawn_local(async move {
                let context = Context::from_slice(&[
                    asset_load_context.as_ref(),
                    config.as_ref(), // TODO: I don't want `Config` in the context long-term
                ]);
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
    let mut asset_load_context = Rc::try_unwrap(asset_load_context)
        .map_err(|_| ())
        .expect("runtime using this context is closed");
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
        self.loader.close(); // TODO: Think about how to actually close the loader, such as cancelling in-progress tasks.
        for task_executor_thread in self.task_executor_threads.drain(..) {
            task_executor_thread.join().unwrap();
        }
        unsafe { queue.drain(&self.gfx.device) };
        unsafe { queue.destroy(&self.gfx.device) };
        unsafe { self.queue_unpark_semaphore.destroy(&self.gfx.device) };
        unsafe {
            Arc::try_unwrap(self.staging.take().unwrap())
                .map_err(|_| ())
                .expect("All threads using staging are now joined")
                .destroy(&self.gfx.device)
        };
    }
}
