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
}

impl AssetLoader {
    pub fn new(gfx: Arc<Base>) -> Self {
        let (cancellation_send, cancellation_receive) = mpsc::channel::<()>();
        let mut queue =
            unsafe { ParallelQueue::new(&gfx.device, gfx.queue_family, gfx.queue, None) };
        let loader = skid_steer::Loader::new();
        let staging = StagingRing::new(&gfx, 32 * 1024 * 1024);

        let join_handle = {
            let gfx = gfx.clone();
            let loader = loader.clone();
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
                        loader.close();
                        unsafe { queue.drive(&gfx.device) };
                        thread::sleep(Duration::from_secs_f32(1.0)); // TODO: We want to park instead, but we need a semaphore for that.
                    }
                });
                unsafe { queue.drain(&gfx.device) };
                unsafe { queue.destroy(&gfx.device) };
                unsafe { staging.destroy(&gfx.device) };
            })
        };

        AssetLoader {
            gfx,
            join_handle: Some(join_handle),
            cancellation_send,
        }
    }
}

impl Drop for AssetLoader {
    fn drop(&mut self) {
        let _ = self.cancellation_send.send(());
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
