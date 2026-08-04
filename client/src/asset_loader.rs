use std::{
    sync::{Arc, mpsc},
    thread::{self, JoinHandle},
};

use ash::vk;
use lahar::{DedicatedMapping, ParallelQueue, TimelineRing};
use skid_steer::Context;

use crate::{Config, graphics::Base};

struct AsyncTimelineRing {
    timeline_ring: TimelineRing,
    size: usize, // TODO: TimelineRing should probably expose this. Otherwise, someone might allocate a buffer of size timeline_ring.capacity() and get undefined behavior.
    max_allocation: usize,
    request_receiver: tokio::sync::mpsc::UnboundedReceiver<AllocationRequest>,
    tick_receiver: tokio::sync::mpsc::UnboundedReceiver<u64>,
    blocked_allocation_request: Option<AllocationRequest>,
}

struct AsyncTimelineRingRequestHandle {
    request_sender: tokio::sync::mpsc::UnboundedSender<AllocationRequest>,
}

struct AsyncTimelineRingTickHandle {
    tick_sender: tokio::sync::mpsc::UnboundedSender<u64>,
}

struct AllocationRequest {
    size: usize,
    align: usize,
    free_at: u64,
    offset_sender: tokio::sync::oneshot::Sender<usize>,
}

impl AsyncTimelineRing {
    fn new(
        size: usize,
    ) -> (
        Self,
        AsyncTimelineRingRequestHandle,
        AsyncTimelineRingTickHandle,
    ) {
        let timeline_ring = TimelineRing::new(size);
        let max_allocation = timeline_ring.capacity();
        let (request_sender, request_receiver) = tokio::sync::mpsc::unbounded_channel();
        let (tick_sender, tick_receiver) = tokio::sync::mpsc::unbounded_channel();
        (
            AsyncTimelineRing {
                timeline_ring,
                size,
                max_allocation,
                request_receiver,
                tick_receiver,
                blocked_allocation_request: None,
            },
            AsyncTimelineRingRequestHandle { request_sender },
            AsyncTimelineRingTickHandle { tick_sender },
        )
    }

    async fn drive(&mut self) {
        loop {
            println!("ParallelQueue loop iteration start");
            tokio::select! {
                biased;
                tick = self.tick_receiver.recv() => {
                    let Some(tick) = tick else {
                        // If the tick sender is closed, that can only mean that we're shutting down.
                        break;
                    };
                    self.apply_tick(tick);
                }
                Some(request) = self.request_receiver.recv(), if self.blocked_allocation_request.is_none() => {
                    self.apply_allocation_request(request);
                }
                else => { break; }
            }
        }
        println!("ParallelQueue loop end");
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

impl AsyncTimelineRingRequestHandle {
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
    timeline_ring: AsyncTimelineRingRequestHandle,
    backing_memory: DedicatedMapping<[u8]>,
    max_allocation: usize,
}

pub struct Allocation<'a> {
    // TODO: Consider encapsulating these fields
    pub bytes: &'a mut [u8],
    pub offset: u64,
}

impl StagingRing {
    fn new(
        gfx: &Base,
        timeline_ring: &AsyncTimelineRing,
        timeline_ring_request_handle: AsyncTimelineRingRequestHandle,
    ) -> Self {
        StagingRing {
            timeline_ring: timeline_ring_request_handle,
            backing_memory: unsafe {
                DedicatedMapping::zeroed_array(
                    &gfx.device,
                    &gfx.memory_properties,
                    vk::BufferUsageFlags::TRANSFER_SRC,
                    timeline_ring.size,
                )
            },
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
}

pub struct ParallelQueueWaiter {
    semaphore: vk::Semaphore,
    queue_unpark_semaphore: HelperSemaphore,
    progress_changed: tokio::sync::Notify,
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

pub struct AssetLoader {
    gfx: Arc<Base>,
    join_handle: Option<JoinHandle<()>>,
    cancellation_send: mpsc::Sender<()>,
    queue_unpark_semaphore: HelperSemaphore,
    loader: skid_steer::Loader,
}

impl AssetLoader {
    pub fn new(gfx: Arc<Base>, config: Arc<Config>) -> Self {
        let (cancellation_send, cancellation_receive) = mpsc::channel::<()>();
        let mut queue =
            unsafe { ParallelQueue::new(&gfx.device, gfx.queue_family, gfx.queue, None) };
        let loader = skid_steer::Loader::new();
        let (mut timeline_ring, timeline_ring_request_handle, timeline_ring_tick_handle) =
            AsyncTimelineRing::new(32 * 1024 * 1024);
        let staging = StagingRing::new(&gfx, &timeline_ring, timeline_ring_request_handle);
        let queue_unpark_semaphore = unsafe { HelperSemaphore::new(&gfx.device) };

        let join_handle = {
            let gfx = gfx.clone();
            let config = config.clone();
            let loader = loader.clone();
            thread::Builder::new()
                .name("asset_loader_driver".to_owned())
                .spawn(move || {
                    let parallel_queue_waiter = ParallelQueueWaiter {
                        semaphore: queue.semaphore(),
                        queue_unpark_semaphore,
                        progress_changed: tokio::sync::Notify::new(),
                    };

                    thread::scope(|s| {
                        let loader = &loader;
                        let staging = &staging;
                        let parallel_queue_waiter = &parallel_queue_waiter;
                        // TODO: Use dynamic number of threads
                        for i in 0..2 {
                            let mut handle = unsafe { queue.handle(&gfx.device) };
                            let gfx: &Base = &gfx;
                            let config: &Config = &config;
                            thread::Builder::new()
                                .name(format!("task_executor_{}", i).to_owned())
                                .spawn_scoped(s, move || {
                                    let runtime = tokio::runtime::LocalRuntime::new().unwrap();
                                    // TODO: This is the wrong pattern for spawning async tasks because concurrency is lost.
                                    runtime.block_on(async {
                                        while let Some(task) = loader.next_task().await {
                                            println!("Found task");
                                            let context = Context::from_slice(&[
                                                gfx,
                                                config, // TODO: I don't want `Config` in the context long-term
                                                &handle,
                                                staging,
                                                parallel_queue_waiter,
                                            ]);
                                            task.run(&context).await;
                                            println!("Asset loaded {}", i);
                                        }
                                        println!("Ending task executor {}", i);
                                        unsafe { handle.destroy(&gfx.device) };
                                    });
                                })
                                .unwrap();
                        }

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
                        loop {
                            // Systems increment the `queue_unpark_semaphore` value when they want to guarantee
                            // that we don't park unless certain things are done. `ParallelQueueWaiter`
                            // wants to ensure that `ParallelQueue::drive` is called, while the cleanup code
                            // wants to ensure `cancellation_recieve` is checked. Therefore, we put these two
                            // operations between reading and waiting on  `queue_unpark_semaphore`
                            let queue_unpark_semaphore_next_value =
                                unsafe { queue_unpark_semaphore.get_next_value(&gfx.device) };
                            unsafe { queue.drive(&gfx.device) };
                            if cancellation_receive.try_recv() != Err(mpsc::TryRecvError::Empty) {
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
                            let _ = timeline_ring_tick_handle.tick_sender.send(semaphore_value);
                        }
                        println!("Cancellation received. Should shut down AssetLoader");
                        loader.close(); // TODO: Think about how to actually close the loader, such as cancelling in-progress tasks.
                        drop(timeline_ring_tick_handle);
                    });
                    println!("All threads shut down");
                    unsafe { queue.drain(&gfx.device) };
                    unsafe { queue.destroy(&gfx.device) };
                    unsafe { queue_unpark_semaphore.destroy(&gfx.device) };
                    unsafe { staging.destroy(&gfx) };
                })
                .unwrap()
        };

        AssetLoader {
            gfx,
            join_handle: Some(join_handle),
            cancellation_send,
            queue_unpark_semaphore,
            loader,
        }
    }

    pub fn load<S: skid_steer::Source>(
        &self,
        source: S,
    ) -> skid_steer::Asset<<S as skid_steer::Source>::Output> {
        self.loader.load(source)
    }
}

impl Drop for AssetLoader {
    fn drop(&mut self) {
        let _ = self.cancellation_send.send(());
        unsafe { self.queue_unpark_semaphore.signal(&self.gfx.device) };
        self.join_handle.take().unwrap().join().unwrap();
    }
}
