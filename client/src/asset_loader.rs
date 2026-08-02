use std::{
    sync::{Arc, Mutex, mpsc},
    thread::{self, JoinHandle},
    time::Duration,
};

use ash::vk;
use lahar::{DedicatedMapping, ParallelQueue, TimelineRing, parallel_queue::Handle};
use skid_steer::Context;

use crate::graphics::Base;

struct AsyncTimelineRing {
    timeline_ring: TimelineRing,
    size: usize, // TODO: TimelineRing should probably expose this. Otherwise, someone might allocate a buffer of size timeline_ring.capacity() and get undefined behavior.
    max_allocation: usize,
    request_sender: tokio::sync::mpsc::UnboundedSender<AllocationRequest>,
    request_receiver: tokio::sync::mpsc::UnboundedReceiver<AllocationRequest>,
    tick_sender: tokio::sync::mpsc::UnboundedSender<u64>,
    tick_receiver: tokio::sync::mpsc::UnboundedReceiver<u64>,
    blocked_allocation_request: Option<AllocationRequest>,
}

struct AsyncTimelineRingHandle {
    request_sender: tokio::sync::mpsc::UnboundedSender<AllocationRequest>,
}

struct AllocationRequest {
    size: usize,
    align: usize,
    free_at: u64,
    offset_sender: tokio::sync::oneshot::Sender<usize>,
}

impl AsyncTimelineRing {
    fn new(size: usize) -> Self {
        let timeline_ring = TimelineRing::new(size);
        let max_allocation = timeline_ring.capacity();
        let (request_sender, request_receiver) = tokio::sync::mpsc::unbounded_channel();
        let (tick_sender, tick_receiver) = tokio::sync::mpsc::unbounded_channel();
        AsyncTimelineRing {
            timeline_ring,
            size,
            max_allocation,
            request_sender,
            request_receiver,
            tick_sender,
            tick_receiver,
            blocked_allocation_request: None,
        }
    }

    fn handle(&self) -> AsyncTimelineRingHandle {
        AsyncTimelineRingHandle {
            request_sender: self.request_sender.clone(),
        }
    }

    async fn drive(&mut self) {
        loop {
            tokio::select! {
                biased;
                Some(tick) = self.tick_receiver.recv() => {
                    self.apply_tick(tick);
                }
                Some(request) = self.request_receiver.recv(), if self.blocked_allocation_request.is_none() => {
                    self.apply_allocation_request(request);
                }
                else => { break; }
            }
        }
    }

    fn apply_tick(&mut self, time: u64) {
        self.timeline_ring.tick(time);

        // Re-attempt allocation request
        if let Some(request) = self.blocked_allocation_request.take() {
            self.apply_allocation_request(request);
        }
    }

    fn apply_allocation_request(&mut self, request: AllocationRequest) {
        if let Some(offset) = self
            .timeline_ring
            .alloc(request.size, request.align, request.free_at)
        {
            let _ = request.offset_sender.send(offset);
        } else {
            self.blocked_allocation_request = Some(request);
        }
    }
}

impl AsyncTimelineRingHandle {
    pub async fn alloc(&self, size: usize, align: usize, free_at: u64) -> usize {
        let (offset_sender, offset_receiver) = tokio::sync::oneshot::channel::<usize>();
        let _ = self.request_sender.send(AllocationRequest {
            size,
            align,
            free_at,
            offset_sender,
        });
        let allocation = offset_receiver.await.expect(
            "all callers of alloc_blocking are dropped before the timeline ring is dropped.",
        );
        allocation // TODO: Does this allocation have a reasonable lifetime?
    }
}

pub struct StagingRing {
    timeline_ring: AsyncTimelineRingHandle,
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
    pub fn new(gfx: &Base, timeline_ring: &AsyncTimelineRing) -> Self {
        StagingRing {
            timeline_ring: timeline_ring.handle(),
            backing_memory: unsafe {
                DedicatedMapping::zeroed_array(
                    &gfx.device,
                    &gfx.memory_properties,
                    vk::BufferUsageFlags::TRANSFER_SRC,
                    timeline_ring.size,
                )
            },
            new_space_available: tokio::sync::Notify::new(),
            max_allocation: timeline_ring.max_allocation.min(isize::MAX as usize),
        }
    }

    /// Safety: The Allocation must be dropped before the staging ring is able to reuse the space taken up by the
    /// allocation. In practice, this means dropping it before calling the `lahar::parallel_queue::Work::end` function.
    pub async unsafe fn alloc_async(
        &self,
        size: usize,
        align: usize,
        free_at: u64,
    ) -> Allocation<'_> {
        assert!(size <= self.max_allocation);
        let offset = self.timeline_ring.alloc(size, align, free_at).await;
        unsafe {
            Allocation {
                bytes: std::slice::from_raw_parts_mut(
                    (self.backing_memory.as_ptr() as *const u8).add(offset) as *mut u8,
                    size,
                ),
                offset: offset.try_into().unwrap(),
            }
        }
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
        let mut timeline_ring = AsyncTimelineRing::new(32 * 1024 * 1024);
        let staging = StagingRing::new(&gfx, &timeline_ring);
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
            thread::Builder::new()
                .name("asset_loader_driver".to_owned())
                .spawn(move || {
                    let parallel_queue_waiter = ParallelQueueWaiter {
                        semaphore: queue.semaphore(),
                        progress_changed: tokio::sync::Notify::new(),
                    };

                    thread::scope(|s| {
                        let loader = &loader;
                        let staging = &staging;
                        let parallel_queue_waiter = &parallel_queue_waiter;
                        // TODO: Use dynamic number of threads
                        for i in 0..2 {
                            let handle = unsafe { queue.handle(&gfx.device) };
                            let gfx: &Base = &gfx;
                            thread::Builder::new()
                                .name(format!("task_executor_{}", i).to_owned())
                                .spawn_scoped(s, move || {
                                    let runtime = tokio::runtime::LocalRuntime::new().unwrap();
                                    runtime.block_on(async {
                                        while let Some(task) = loader.next_task().await {
                                            let mut context = Context::new();
                                            context.insert::<Base>(gfx);
                                            context.insert::<Handle>(&handle);
                                            context.insert::<StagingRing>(staging);
                                            context.insert::<ParallelQueueWaiter>(
                                                parallel_queue_waiter,
                                            );
                                            task.run(&context).await;
                                        }
                                    });
                                })
                                .unwrap();
                        }

                        let tick_sender = timeline_ring.tick_sender.clone();
                        thread::Builder::new()
                            .name("timeline_ring_driver".to_owned())
                            .spawn_scoped(s, move || {
                                let runtime = tokio::runtime::LocalRuntime::new().unwrap();
                                runtime.block_on(async {
                                    timeline_ring.drive().await;
                                })
                            })
                            .unwrap();

                        // Driver thread
                        while cancellation_receive.try_recv() == Err(mpsc::TryRecvError::Empty) {
                            unsafe { queue.drive(&gfx.device) };
                            let semaphore_value =
                                unsafe { queue.park(&gfx.device, cancellation_semaphore, 1) };
                            parallel_queue_waiter.progress_changed.notify_waiters();
                            let _ = tick_sender.send(semaphore_value);
                        }
                    });
                    loader.close(); // TODO: Think about how to actually close the loader, such as cancelling in-progress tasks.
                    unsafe { queue.drain(&gfx.device) };
                    unsafe { queue.destroy(&gfx.device) };
                    unsafe { gfx.device.destroy_semaphore(cancellation_semaphore, None) };
                    unsafe { staging.destroy(&gfx) };
                })
                .unwrap()
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
