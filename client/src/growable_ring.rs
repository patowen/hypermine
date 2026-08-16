use std::{collections::VecDeque, mem, ptr::NonNull, sync::Mutex};

use ash::{Device, vk};
use lahar::{DedicatedBuffer, TimelineRing};

use crate::graphics::Base;

pub struct GrowableRing {
    state: Mutex<State>,
    old: Mutex<VecDeque<(u64, DedicatedBuffer)>>,
}

pub struct Allocation<T> {
    pub buffer: vk::Buffer,
    pub offset: u64,
    pub size: u64,
    pub pointer: NonNull<T>,
}

impl GrowableRing {
    pub fn new(gfx: &Base, capacity: usize) -> Self {
        Self {
            state: Mutex::new(State::new(gfx, capacity)),
            old: Mutex::new(VecDeque::new()),
        }
    }

    pub unsafe fn destroy(&mut self, device: &Device) {
        unsafe {
            self.state.get_mut().unwrap().memory.destroy(device);
            for (_, mut buffer) in self.old.get_mut().unwrap().drain(..) {
                buffer.destroy(device);
            }
        }
    }

    /// The returned pointer must not be dereffed after `tick` is called with `free_at`
    pub fn alloc<T>(&self, gfx: &Base, count: usize, align: usize, free_at: u64) -> Allocation<T> {
        let size = count * mem::size_of::<T>();
        let (buffer, offset, mapping) = {
            let mut state = self.state.lock().unwrap();
            let offset = if let Some(offset) =
                state
                    .alloc
                    .alloc(size, align.max(mem::align_of::<T>()), free_at)
            {
                offset
            } else {
                let old = mem::replace(&mut *state, State::new(gfx, size * 2));
                // free_at is
                self.old
                    .lock()
                    .unwrap()
                    .push_back((old.free_at, old.memory));
                state
                    .alloc
                    .alloc(size, align.max(mem::align_of::<T>()), free_at)
                    .expect("alloc failed after growing")
            };
            state.free_at = state.free_at.max(free_at);
            (state.memory.handle, offset, unsafe {
                NonNull::new_unchecked(state.mapping.as_ptr().add(offset).cast())
            })
        };
        Allocation {
            buffer,
            offset: offset as u64,
            size: size as u64,
            pointer: mapping,
        }
    }

    pub unsafe fn tick(&self, device: &Device, time: u64) {
        unsafe {
            self.state.lock().unwrap().alloc.tick(time);
            let mut old = self.old.lock().unwrap();
            while let Some(&(t, _)) = old.front() {
                if t > time {
                    break;
                }
                let (_, mut buffer) = old.pop_front().unwrap();
                buffer.destroy(device);
            }
        }
    }
}

struct State {
    alloc: TimelineRing,
    memory: DedicatedBuffer,
    mapping: NonNull<u8>,
    /// Largest timeline value any allocation may be in use before
    free_at: u64,
}

impl State {
    fn new(gfx: &Base, capacity: usize) -> Self {
        let alloc = TimelineRing::new(capacity);
        let memory = unsafe {
            DedicatedBuffer::new(
                &gfx.device,
                &gfx.memory_properties,
                &vk::BufferCreateInfo::default()
                    .size(alloc.capacity() as vk::DeviceSize + 1)
                    .usage(vk::BufferUsageFlags::TRANSFER_SRC)
                    .sharing_mode(vk::SharingMode::EXCLUSIVE),
                vk::MemoryPropertyFlags::HOST_VISIBLE | vk::MemoryPropertyFlags::HOST_COHERENT,
            )
        };
        let mapping = unsafe {
            gfx.set_name(memory.handle, c"staging");
            NonNull::new_unchecked(
                gfx.device
                    .map_memory(
                        memory.memory,
                        0,
                        vk::WHOLE_SIZE,
                        vk::MemoryMapFlags::default(),
                    )
                    .unwrap()
                    .cast(),
            )
        };
        Self {
            alloc,
            memory,
            mapping,
            free_at: 0,
        }
    }
}

unsafe impl Send for GrowableRing {}
unsafe impl Sync for GrowableRing {}
