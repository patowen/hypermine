use std::{
    sync::{Arc, Mutex},
    thread,
};

use ash::vk;
use lahar::{DedicatedMapping, ParallelQueue, TimelineRing};
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
    pub fn new(size: usize) {
        StagingRing {
            timeline_ring: Mutex::new(TimelineRing::new(size)),
            backing_memory: unsafe {
                DedicatedMapping::zeroed_array(todo!(), todo!(), todo!(), size)
            },
        };
    }

    /// Safety: The allocation is not allowed to be freed before the returned reference's lifetime ends. (TODO: Explain this better)
    pub unsafe fn alloc(&self, size: usize, align: usize, free_at: u64) -> Option<Allocation> {
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
}

fn commence(gfx: Arc<Base>) {
    let mut queue = unsafe { ParallelQueue::new(&gfx.device, gfx.queue_family, gfx.queue, None) };
    let loader = skid_steer::Loader::new();

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
        s.spawn(move || {
            unsafe { queue.drive(&gfx.device) };
        });
    });
}
