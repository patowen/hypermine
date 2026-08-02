use std::{
    sync::{Arc, Mutex, mpsc},
    thread::{self, JoinHandle},
    time::Duration,
};

use ash::vk;
use lahar::{DedicatedMapping, ParallelQueue, TimelineRing, parallel_queue::Handle};
use skid_steer::Context;

use crate::graphics::Base;

pub struct StagingRing {
    timeline_ring: Mutex<TimelineRing>,
    backing_memory: DedicatedMapping<[u8]>,
    new_space_available: tokio::sync::Notify,
    max_allocation: usize,
}

pub struct Allocation<'a> {
    // TODO: Consider encapsulating these fields
    pub bytes: &'a mut [u8],
    pub offset: u64,
}

impl StagingRing {
    pub fn new(gfx: &Base, size: usize) -> Self {
        let timeline_ring = TimelineRing::new(size);
        let max_allocation = timeline_ring.capacity().min(isize::MAX as usize);
        StagingRing {
            timeline_ring: Mutex::new(timeline_ring),
            backing_memory: unsafe {
                DedicatedMapping::zeroed_array(
                    &gfx.device,
                    &gfx.memory_properties,
                    vk::BufferUsageFlags::TRANSFER_SRC,
                    size,
                )
            },
            new_space_available: tokio::sync::Notify::new(),
            max_allocation,
        }
    }

    pub async unsafe fn alloc_blocking(
        &self,
        size: usize,
        align: usize,
        free_at: u64,
    ) -> Allocation<'_> {
        // TODO: Instead of blocking in this haphazard way, try growing the timeline ring instead if it's too small.
        // If we decide blocking is better, try queueing allocations instead of racing for a mutex for better fairness.
        assert!(size <= self.max_allocation);
        loop {
            let notified = self.new_space_available.notified();
            if let Some(alloc) = unsafe { self.try_alloc(size, align, free_at) } {
                return alloc;
            }
            notified.await;
        }
    }

    /// Safety: The Allocation must be dropped before the staging ring is able to reuse the space taken up by the
    /// allocation. In practice, this means dropping it before calling the `lahar::parallel_queue::Work::end` function.
    pub unsafe fn try_alloc(
        &self,
        size: usize,
        align: usize,
        free_at: u64,
    ) -> Option<Allocation<'_>> {
        // TODO: If we decide to queue allocations, then the raw allocation function will no longer need a mutex.
        assert!(size <= self.max_allocation);
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

    /// Safety: gfx needs to be the same as the Base passed to `new`
    pub unsafe fn destroy(mut self, gfx: &Base) {
        unsafe { self.backing_memory.destroy(&gfx.device) };
    }

    // TODO: tick
}

pub struct ParallelQueueWaiter {
    semaphore: vk::Semaphore,
    progress_changed: tokio::sync::Notify,
}

impl ParallelQueueWaiter {
    pub async fn wait_for_semaphore(&self, device: &ash::Device, value: u64) {
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

struct AssetLoader {
    gfx: Arc<Base>,
    join_handle: Option<JoinHandle<()>>,
    cancellation_send: mpsc::Sender<()>,
    cancellation_semaphore: vk::Semaphore,
}

impl AssetLoader {
    pub fn new(gfx: Arc<Base>) -> Self {
        let (cancellation_send, cancellation_receive) = mpsc::channel::<()>();
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
            thread::spawn(move || {
                let parallel_queue_waiter = ParallelQueueWaiter {
                    semaphore: queue.semaphore(),
                    progress_changed: tokio::sync::Notify::new(),
                };

                thread::scope(|s| {
                    let loader = &loader;
                    let staging = &staging;
                    let parallel_queue_waiter = &parallel_queue_waiter;
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
                                    context.insert::<ParallelQueueWaiter>(parallel_queue_waiter);
                                    task.run(&context).await;
                                }
                            });
                        });
                    }

                    // Driver thread
                    while cancellation_receive.try_recv() == Err(mpsc::TryRecvError::Empty) {
                        unsafe { queue.drive(&gfx.device) };
                        unsafe { queue.park(&gfx.device, cancellation_semaphore, 1) };
                        parallel_queue_waiter.progress_changed.notify_waiters();
                    }
                });
                loader.close(); // TODO: Think about how to actually close the loader, such as cancelling in-progress tasks.
                unsafe { queue.drain(&gfx.device) };
                unsafe { queue.destroy(&gfx.device) };
                unsafe { gfx.device.destroy_semaphore(cancellation_semaphore, None) };
                unsafe { staging.destroy(&gfx) };
            })
        };

        AssetLoader {
            gfx,
            join_handle: Some(join_handle),
            cancellation_send,
            cancellation_semaphore,
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

// TODO: This was an earlier draft. I should delete this.
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
    unsafe { staging.destroy(&gfx) };
}
