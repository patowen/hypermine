use std::{
    sync::{Arc, Mutex, atomic::AtomicU64, mpsc},
    thread::{self, JoinHandle},
    time::Duration,
};

use ash::vk::{self, SemaphoreSignalInfo};
use lahar::{DedicatedMapping, ParallelQueue, TimelineRing, parallel_queue::Handle};
use skid_steer::Context;

use crate::graphics::Base;

pub struct StagingRing {
    timeline_ring: Mutex<TimelineRing>,
    backing_memory: DedicatedMapping<[u8]>,
}

pub struct Allocation<'a> {
    // TODO: Consider encapsulating these fields
    pub bytes: &'a mut [u8],
    pub offset: u64,
}

impl StagingRing {
    pub fn new(gfx: &Arc<Base>, size: usize) -> Self {
        StagingRing {
            timeline_ring: Mutex::new(TimelineRing::new(size)),
            backing_memory: unsafe {
                DedicatedMapping::zeroed_array(
                    &gfx.device,
                    &gfx.memory_properties,
                    vk::BufferUsageFlags::TRANSFER_SRC,
                    size,
                )
            },
        }
    }

    /// Safety: The allocation is not allowed to be freed before the returned reference's lifetime ends. (TODO: Explain this better)
    pub unsafe fn alloc(&self, size: usize, align: usize, free_at: u64) -> Option<Allocation<'_>> {
        // TODO: Need a way to wait for an allocation
        let offset = self
            .timeline_ring
            .lock()
            .unwrap()
            .alloc(size, align, free_at)?;
        Some(unsafe {
            Allocation {
                bytes: std::slice::from_raw_parts_mut(
                    (self.backing_memory.as_ptr() as *const u8).add(offset) as *mut u8,
                    size,
                ),
                offset: offset.try_into().unwrap(),
            }
        })
    }

    pub fn buffer(&self) -> vk::Buffer {
        self.backing_memory.buffer()
    }

    /// Safety: Device needs to be the same as the device passed to `new`
    pub unsafe fn destroy(mut self, device: &ash::Device) {
        unsafe { self.backing_memory.destroy(device) };
    }
}

struct AssetLoader {
    gfx: Arc<Base>,
    join_handle: Option<JoinHandle<()>>,
    cancellation_send: mpsc::Sender<()>,
    cancellation_semaphore: vk::Semaphore,
    progress_changed: Arc<tokio::sync::Notify>,
    current_progress: Arc<AtomicU64>,
}

impl AssetLoader {
    pub fn new(gfx: Arc<Base>) -> Self {
        let (cancellation_send, cancellation_receive) = mpsc::channel::<()>();
        let progress_changed = Arc::new(tokio::sync::Notify::new());
        let current_progress = Arc::new(AtomicU64::new(0));
        let mut queue =
            unsafe { ParallelQueue::new(&gfx.device, gfx.queue_family, gfx.queue, None) };
        let loader = skid_steer::Loader::new();
        let staging = StagingRing::new(&gfx, 32 * 1024 * 1024);
        let cancellation_semaphore = unsafe {
            gfx.device.create_semaphore(
                &vk::SemaphoreCreateInfo::default().push_next(
                    &mut vk::SemaphoreTypeCreateInfo::default()
                        .semaphore_type(vk::SemaphoreType::TIMELINE)
                        .initial_value(0),
                ),
                None,
            )
        }
        .unwrap();

        let join_handle = {
            let gfx = gfx.clone();
            let loader = loader.clone();
            let current_progress = current_progress.clone();
            let progress_changed = progress_changed.clone();
            thread::spawn(move || {
                thread::scope(|s| {
                    let loader = &loader;
                    let staging = &staging;
                    // TODO: Use dynamic number of threads
                    for _ in 0..2 {
                        let handle = unsafe { queue.handle(&gfx.device) };
                        let gfx: &Base = &gfx;
                        s.spawn(move || {
                            let runtime = tokio::runtime::LocalRuntime::new().unwrap();
                            runtime.block_on(async {
                                while let Some(task) = loader.next_task().await {
                                    let mut context = Context::new();
                                    context.insert::<Base>(gfx);
                                    context.insert::<Handle>(&handle);
                                    context.insert::<StagingRing>(staging);
                                    task.run(&context).await;
                                }
                            });
                        });
                    }

                    // Driver thread
                    while cancellation_receive.try_recv() == Err(mpsc::TryRecvError::Empty) {
                        unsafe { queue.drive(&gfx.device) };
                        let timeline_value =
                            unsafe { queue.park(&gfx.device, cancellation_semaphore, 1) };
                        current_progress
                            .store(timeline_value, std::sync::atomic::Ordering::Relaxed);
                        progress_changed.notify_waiters();
                    }
                });
                loader.close(); // TODO: Think about how to actually close the loader, such as cancelling in-progress tasks.
                unsafe { queue.drain(&gfx.device) };
                unsafe { queue.destroy(&gfx.device) };
                unsafe { gfx.device.destroy_semaphore(cancellation_semaphore, None) };
                unsafe { staging.destroy(&gfx.device) };
            })
        };

        AssetLoader {
            gfx,
            join_handle: Some(join_handle),
            cancellation_send,
            cancellation_semaphore,
            progress_changed,
            current_progress: Arc::new(AtomicU64::new(0)),
        }
    }

    pub async fn wait_for_parallel_queue_work_completion(&self, value: u64) {
        let mut notified = self.progress_changed.notified();
        loop {
            if self
                .current_progress
                .load(std::sync::atomic::Ordering::Relaxed)
                > value /* TODO: Check for off-by-one error. What do we actually want to wait on? */
            {
                return;
            }
            /*
                Justification for why the following await shouldn't deadlock:
                We just need to wait until `notify_waiters()` is called since `notified` was set, so the main question is:
                    Is there any scenario in which the final call to `notify_waiters()` was already missed, but the atomic, when loaded, was not at the right value for exiting the loop?
                Initializing and awaiting `notified` are both acquire operations on the state of the `Notify`, while `notify_waiters` is a release operation on the same state.
                The only way `notified` can miss the `notify_waiters()` call is if it saw the state that `notify_waiters()` had already produced, meaning that `notify_waiters()` happened_before the setting of `notified`.
                However, this would also mean that `current_progress` being set to the loop-exiting value must have happened_before `current_progress` was read, making this deadlock scenario impossible.

                TODO: Is there a simpler mental model for this? Gaining confidence in the correct usage of atomics seems incredibly difficult.
             */
            notified.await;
            notified = self.progress_changed.notified();
        }
    }
}

impl Drop for AssetLoader {
    fn drop(&mut self) {
        let _ = self.cancellation_send.send(());
        unsafe {
            self.gfx
                .device
                .signal_semaphore(
                    &vk::SemaphoreSignalInfo::default()
                        .semaphore(self.cancellation_semaphore)
                        .value(1),
                )
                .unwrap()
        };
        self.join_handle.take().unwrap().join().unwrap();
    }
}

fn commence(gfx: Arc<Base>) {
    let mut queue = unsafe { ParallelQueue::new(&gfx.device, gfx.queue_family, gfx.queue, None) };
    let loader = skid_steer::Loader::new();
    let staging = StagingRing::new(&gfx, 32 * 1024 * 1024);

    thread::scope(|s| {
        // TODO: Use dynamic number of threads
        for _ in 0..2 {
            let handle = unsafe { queue.handle(&gfx.device) };
            let gfx = gfx.clone();
            let loader = loader.clone();
            s.spawn(move || {
                let runtime = tokio::runtime::LocalRuntime::new().unwrap();
                runtime.block_on(async {
                    while let Some(task) = loader.next_task().await {
                        let mut context = Context::new();
                        context.insert(&gfx);
                        context.insert(&handle);
                        task.run(&context).await;
                    }
                });
            });
        }

        // Driver thread
        let gfx = gfx.clone();
        loop {
            unsafe { queue.drive(&gfx.device) };
            thread::sleep(Duration::from_secs_f32(1.0));
        }
    });

    unsafe { queue.destroy(&gfx.device) };
    unsafe { staging.destroy(&gfx.device) };
}
