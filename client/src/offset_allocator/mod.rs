//! This module is a modified version of the https://github.com/pcwalton/offset-allocator crate,
//! updated for improved readability and optimized for Hypermine's use case.
//! It is based on commit 9449d4d81be91396938b05894fad443f63093077.

/* The license notice of the original offset-allocator crate this came from is as follows:

MIT License

Copyright (c) 2023 Sebastian Aaltonen, Patrick Walton

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

mod bins_map;
mod node_index;
mod small_float;

use std::fmt::{Debug, Formatter};

use slab::Slab;
use tracing::debug;

use crate::offset_allocator::{
    bins_map::BinsMap,
    node_index::{NodeIndex, NodeIndexOption},
    small_float::{SmallFloat, SmallFloatMap},
};

/// An allocator that manages a single contiguous chunk of space and hands out
/// portions of it as requested. Since this allocator does not support alignment, it is recommended
/// to interpret these allocations in whatever unit is most convenient, which will likely not be "bytes"
pub struct Allocator<NI: NodeIndex = u32> {
    /// The total size of the buffer
    size: u32,
    /// The maximum number of "nodes", or continuous blocks the allocator can handle. The actual supported number of allocations is less than this.
    max_nodes: u32,
    /// The total amount of remaining available space in the buffer. Fragmentation and rounding means that an allocation of this size is not always possible,
    /// but as long as this is non-zero, and `max_nodes` isn't exceeded, it's always possible to create an allocation of size 1.
    free_storage: u32,
    /// A [`BinsMap`] that keeps track of all nodes that are not part of an existing allocation
    bins_map: BinsMap<NI>,
    /// Maintains the mapping from [`NodeIndex`] to [`Node`]
    nodes: NodeSlab<NI>,
}

/// A single allocation.
#[derive(Clone, Copy)]
pub struct Allocation<NI: NodeIndex = u32> {
    /// The location of this allocation within the buffer.
    pub offset: u32,
    /// The node index associated with this allocation.
    metadata: NI,
}

/// Provides a summary of the state of the allocator, including space remaining.
#[derive(Debug)]
pub struct StorageReport {
    /// The amount of free space left.
    pub total_free_space: u32,
    /// The maximum potential size of a single contiguous allocation.
    pub largest_free_region: u32,
}

/// Provides a detailed accounting of each bin within the allocator.
#[derive(Debug, Default)]
pub struct StorageReportFull {
    /// Each bin within the allocator.
    pub free_regions: SmallFloatMap<StorageReportFullRegion>,
}

/// A detailed accounting of each allocator bin.
#[derive(Clone, Copy, Debug, Default)]
pub struct StorageReportFullRegion {
    /// The size of the bin, in units.
    pub size: u32,
    /// The number of allocations in the bin.
    pub count: u32,
}

#[derive(Clone, Copy)]
struct Node<NI: NodeIndex = u32> {
    /// The offset of the node in the buffer
    data_offset: u32,
    /// The size of the node in the buffer
    data_size: u32,
    /// Nodes representing free space are added to bins based on their size. Each bin can store an arbitrary number of nodes,
    /// so we used a linked list. This stores the previous node in the bin. This field is meaningless when the node is used in an active allocation.
    bin_list_prev: NodeIndexOption<NI>,
    /// Nodes representing free space are added to bins based on their size. Each bin can store an arbitrary number of nodes,
    /// so we used a linked list. This stores the next node in the bin. This field is meaningless when the node is used in an active allocation.
    bin_list_next: NodeIndexOption<NI>,
    /// The entire buffer is split up into several nodes, some marking an allocation and others marking free space.
    /// Neighboring nodes in this buffer point to each other in a linked list. This field stores the index of the previous neighboring node.
    neighbor_prev: NodeIndexOption<NI>,
    /// The entire buffer is split up into several nodes, some marking an allocation and others marking free space.
    /// Neighboring nodes in this buffer point to each other in a linked list. This field stores the index of the next neighboring node.
    neighbor_next: NodeIndexOption<NI>,
    /// Whether the node is used in an active allocation
    used: bool, // Note: One possible enhancement to reduce the size of `Node` is to merge this with another field as a bit flag.
}

/// A slab that stores [`Node`]s by their [`NodeIndex`]. The maximum capacity of this slab depends on the size of the [`NodeIndex`].
struct NodeSlab<NI: NodeIndex>(Slab<Node<NI>>);

impl<NI: NodeIndex> NodeSlab<NI> {
    /// Construct a new, empty `NodeSlab`
    #[inline]
    pub fn new() -> Self {
        // TODO: To avoid potential pauses, we may want to pre-allocate the maximum number of allocations we need.
        NodeSlab(Slab::new())
    }

    /// Return the number of stored nodes
    #[inline]
    pub fn len(&self) -> u32 {
        self.0.len() as u32
    }

    /// Insert a node into the slab, returning the index associated with it
    #[inline]
    pub fn insert(&mut self, node: Node<NI>) -> NI {
        assert!(self.len() != u32::MAX);
        NI::from_u32(self.0.insert(node) as u32)
    }

    /// Remove and return the node associated with the index
    #[inline]
    pub fn remove(&mut self, index: NI) -> Node<NI> {
        self.0.remove(index.to_usize())
    }
}

impl<NI: NodeIndex> std::ops::Index<NI> for NodeSlab<NI> {
    type Output = Node<NI>;

    #[inline]
    fn index(&self, index: NI) -> &Self::Output {
        &self.0[index.to_usize()]
    }
}

impl<NI: NodeIndex> std::ops::IndexMut<NI> for NodeSlab<NI> {
    #[inline]
    fn index_mut(&mut self, index: NI) -> &mut Self::Output {
        &mut self.0[index.to_usize()]
    }
}

impl<NI: NodeIndex> Allocator<NI> {
    /// Creates a new allocator, managing a contiguous block of memory of `size`
    /// units, with the maximum allocations set as high as possible.
    pub fn new(size: u32) -> Self {
        Allocator::with_max_nodes(size, NI::NUM_VALID)
    }

    /// Creates a new allocator, managing a contiguous block of memory of `size`
    /// units, with the given number of maximum nodes.
    /// 
    /// Note that even if no memory is freed, the maximum number of allocations
    /// allowed is 1 less than the maximum number of nodes, since a node is needed
    /// to keep track of the remaining free space. If memory is freed, due to fragmentation,
    /// it is not guaranteed that another allocation will become available.
    ///
    /// Note also that the maximum number of nodes must be at most
    /// [`NodeIndex::NUM_VALID`] at at least 1. If this restriction is violated, this
    /// constructor will panic.
    pub fn with_max_nodes(size: u32, max_nodes: u32) -> Self {
        assert!(max_nodes > 0);
        assert!(max_nodes <= NI::NUM_VALID);

        let mut this = Self {
            size,
            max_nodes,
            free_storage: 0,
            bins_map: BinsMap::default(),
            nodes: NodeSlab::new(),
        };
        this.insert_node_into_bin(size, 0);
        this
    }

    /// Clears out all allocations.
    pub fn reset(&mut self) {
        *self = Self::with_max_nodes(self.size, self.max_nodes);
    }

    /// Allocates a block of `size` elements and returns its allocation.
    ///
    /// If there's not enough contiguous space for this allocation, returns
    /// None.
    pub fn allocate(&mut self, size: u32) -> Option<Allocation<NI>> {
        // Out of allocations?
        if self.nodes.len() >= self.max_nodes {
            return None;
        }

        // Round up when finding the bin index to ensure that any node in that bin can hold the allocation
        let min_bin_index = SmallFloat::from_u32_round_up(size);
        let bin_index = self.bins_map.min_occupied_since(min_bin_index)?;

        // Pop the top node of the bin from the linked list
        let node_index = self.bins_map[bin_index].unwrap();
        let node = &mut self.nodes[node_index];
        let node_total_size = node.data_size;
        node.data_size = size;
        node.used = true;
        self.bins_map
            .replace_bin_node(bin_index, node.bin_list_next);
        if let Some(bin_list_next) = node.bin_list_next.to_option() {
            self.nodes[bin_list_next].bin_list_prev = NodeIndexOption::NONE;
        }
        self.free_storage -= node_total_size;
        debug!(
            "Free storage: {} (-{}) (allocate)",
            self.free_storage, node_total_size
        );

        // Push back remainder N elements to a (usually) lower bin
        let remainder_size = node_total_size - size;
        if remainder_size > 0 {
            let Node {
                data_offset,
                neighbor_next,
                ..
            } = self.nodes[node_index];

            let new_node_index = self.insert_node_into_bin(remainder_size, data_offset + size);

            // Link nodes next to each other so that we can merge them later if both are free
            // And update the old next neighbor to point to the new node (in middle)
            let node = &mut self.nodes[node_index];
            if let Some(neighbor_next) = node.neighbor_next.to_option() {
                self.nodes[neighbor_next].neighbor_prev = NodeIndexOption::some(new_node_index);
            }
            self.nodes[new_node_index].neighbor_prev = NodeIndexOption::some(node_index);
            self.nodes[new_node_index].neighbor_next = neighbor_next;
            self.nodes[node_index].neighbor_next = NodeIndexOption::some(new_node_index);
        }

        let node = &mut self.nodes[node_index];
        Some(Allocation {
            offset: node.data_offset,
            metadata: node_index,
        })
    }

    /// Frees an allocation, returning the data to the heap.
    ///
    /// If the allocation has already been freed, the behavior is unspecified.
    /// It may or may not panic. Note that the memory safety of the allocator *itself* will be
    /// uncompromised, even on double free.
    pub fn free(&mut self, allocation: Allocation<NI>) {
        let node_index = allocation.metadata;

        // Merge with neighbors…
        let Node {
            data_offset: mut offset,
            data_size: mut size,
            used,
            ..
        } = self.nodes[node_index];

        // Double delete check
        assert!(used);

        if let Some(neighbor_prev) = self.nodes[node_index].neighbor_prev.to_option() {
            if !self.nodes[neighbor_prev].used {
                // Previous (contiguous) free node: Change offset to previous
                // node offset. Sum sizes
                let prev_node = &self.nodes[neighbor_prev];
                offset = prev_node.data_offset;
                size += prev_node.data_size;

                let prev_node = &self.nodes[neighbor_prev];
                debug_assert_eq!(prev_node.neighbor_next, NodeIndexOption::some(node_index));
                self.nodes[node_index].neighbor_prev = prev_node.neighbor_prev;
                self.remove_node_from_bin(neighbor_prev);
            }
        }

        if let Some(neighbor_next) = self.nodes[node_index].neighbor_next.to_option() {
            if !self.nodes[neighbor_next].used {
                // Next (contiguous) free node: Offset remains the same. Sum
                // sizes.
                let next_node = &self.nodes[neighbor_next];
                size += next_node.data_size;

                let next_node = &self.nodes[neighbor_next];
                debug_assert_eq!(next_node.neighbor_prev, NodeIndexOption::some(node_index));
                self.nodes[node_index].neighbor_next = next_node.neighbor_next;
                self.remove_node_from_bin(neighbor_next);
            }
        }

        let Node {
            neighbor_next,
            neighbor_prev,
            ..
        } = self.nodes[node_index];

        self.nodes.remove(node_index);

        // Insert the (combined) free node to bin
        let combined_node_index = self.insert_node_into_bin(size, offset);

        // Connect neighbors with the new combined node
        if let Some(neighbor_next) = neighbor_next.to_option() {
            self.nodes[combined_node_index].neighbor_next = NodeIndexOption::some(neighbor_next);
            self.nodes[neighbor_next].neighbor_prev = NodeIndexOption::some(combined_node_index);
        }
        if let Some(neighbor_prev) = neighbor_prev.to_option() {
            self.nodes[combined_node_index].neighbor_prev = NodeIndexOption::some(neighbor_prev);
            self.nodes[neighbor_prev].neighbor_next = NodeIndexOption::some(combined_node_index);
        }
    }

    /// Creates a new free [`Node`] and inserts it at the head of the appropriate bin. Note that the caller of this
    /// function is responsible for linking node in the "neighbor" linked list.
    fn insert_node_into_bin(&mut self, size: u32, data_offset: u32) -> NI {
        // Round down when finding the bin index to ensure that the node being put in that bin can hold any allocation associated with that bin
        let bin_index = SmallFloat::from_u32_round_down(size);

        // Create a new node and insert on top of the bin linked list
        let top_node_index = self.bins_map[bin_index];
        let node_index = self.nodes.insert(Node {
            data_offset,
            data_size: size,
            bin_list_prev: NodeIndexOption::NONE,
            bin_list_next: top_node_index,
            neighbor_prev: NodeIndexOption::NONE,
            neighbor_next: NodeIndexOption::NONE,
            used: false,
        });
        if let Some(top_node_index) = top_node_index.to_option() {
            self.nodes[top_node_index].bin_list_prev = NodeIndexOption::some(node_index);
        }
        self.bins_map
            .replace_bin_node(bin_index, NodeIndexOption::some(node_index));

        self.free_storage += size;
        debug!(
            "Free storage: {} (+{}) (insert_node_into_bin)",
            self.free_storage, size
        );
        node_index
    }

    /// Deletes a [`Node`], removing it from the bin. Note that the caller of this
    /// function is responsible for fixing up links in the "neighbor" linked list, and it is recommended
    /// that this fixup occur before this function is called.
    fn remove_node_from_bin(&mut self, node_index: NI) {
        // Copy the node to work around borrow check.
        let node = self.nodes[node_index];

        match node.bin_list_prev.to_option() {
            Some(bin_list_prev) => {
                // Easy case: We have previous node. Just remove this node from the middle of the list.
                self.nodes[bin_list_prev].bin_list_next = node.bin_list_next;
                if let Some(bin_list_next) = node.bin_list_next.to_option() {
                    self.nodes[bin_list_next].bin_list_prev = node.bin_list_prev;
                }
            }
            None => {
                // Hard case: We are the first node in a bin. Find the bin.

                // Round down when finding the bin index to ensure consistency with `insert_node_into_bin`
                let bin_index = SmallFloat::from_u32_round_down(node.data_size);

                self.bins_map
                    .replace_bin_node(bin_index, node.bin_list_next);
                if let Some(bin_list_next) = node.bin_list_next.to_option() {
                    self.nodes[bin_list_next].bin_list_prev = NodeIndexOption::NONE;
                }
            }
        }

        self.nodes.remove(node_index);

        self.free_storage -= node.data_size;
        debug!(
            "Free storage: {} (-{}) (remove_node_from_bin)",
            self.free_storage, node.data_size
        );
    }

    /// Returns the *used* size of an allocation.
    ///
    /// For this allocator, this always equals the size requested at allocation time.
    pub fn allocation_size(&self, allocation: Allocation<NI>) -> u32 {
        self.nodes[allocation.metadata].data_size
    }

    /// Returns a structure containing the amount of free space remaining, as
    /// well as the largest amount that can be allocated at once.
    pub fn storage_report(&self) -> StorageReport {
        if self.nodes.len() >= self.max_nodes {
            // Out of allocations? -> Zero free space
            return StorageReport {
                total_free_space: 0,
                largest_free_region: 0,
            };
        }

        let largest_free_region = self.bins_map.max_occupied().map_or(0, |x| x.to_u32());
        debug_assert!(self.free_storage >= largest_free_region);

        StorageReport {
            total_free_space: self.free_storage,
            largest_free_region,
        }
    }

    /// Returns detailed information about the number of allocations in each bin.
    pub fn storage_report_full(&self) -> StorageReportFull {
        let mut report = StorageReportFull::default();
        for i in SmallFloat::values() {
            let mut count = 0;
            let mut maybe_node_index = self.bins_map[i];
            while let Some(node_index) = maybe_node_index.to_option() {
                maybe_node_index = self.nodes[node_index].bin_list_next;
                count += 1;
            }
            report.free_regions[i] = StorageReportFullRegion {
                size: i.to_u32(),
                count,
            }
        }
        report
    }
}

impl<NI: NodeIndex> Debug for Allocator<NI> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        self.storage_report().fmt(f)
    }
}

/// Returns the minimum allocator size needed to hold an object of the given size.
pub fn min_allocator_size(needed_object_size: u32) -> u32 {
    SmallFloat::from_u32_round_up(needed_object_size).to_u32()
}

#[cfg(test)]
mod tests {
    use crate::offset_allocator;

    use super::*;

    #[test]
    fn basic_offset_allocator() {
        let mut allocator: Allocator = Allocator::new(1024 * 1024 * 256);
        let a = allocator.allocate(1337).unwrap();
        let offset = a.offset;
        assert_eq!(offset, 0);
        allocator.free(a);
    }

    #[test]
    fn allocate_offset_allocator_simple() {
        let mut allocator: Allocator<u32> = Allocator::new(1024 * 1024 * 256);

        // Free merges neighbor empty nodes. Next allocation should also have offset = 0
        let a = allocator.allocate(0).unwrap();
        assert_eq!(a.offset, 0);

        let b = allocator.allocate(1).unwrap();
        assert_eq!(b.offset, 0);

        let c = allocator.allocate(123).unwrap();
        assert_eq!(c.offset, 1);

        let d = allocator.allocate(1234).unwrap();
        assert_eq!(d.offset, 124);

        allocator.free(a);
        allocator.free(b);
        allocator.free(c);
        allocator.free(d);

        // End: Validate that allocator has no fragmentation left. Should be 100% clean.
        let validate_all = allocator.allocate(1024 * 1024 * 256).unwrap();
        assert_eq!(validate_all.offset, 0);
        allocator.free(validate_all);
    }

    #[test]
    fn allocate_offset_allocator_merge_trivial() {
        let mut allocator: Allocator<u32> = Allocator::new(1024 * 1024 * 256);

        // Free merges neighbor empty nodes. Next allocation should also have offset = 0
        let a = allocator.allocate(1337).unwrap();
        assert_eq!(a.offset, 0);
        allocator.free(a);

        let b = allocator.allocate(1337).unwrap();
        assert_eq!(b.offset, 0);
        allocator.free(b);

        // End: Validate that allocator has no fragmentation left. Should be 100% clean.
        let validate_all = allocator.allocate(1024 * 1024 * 256).unwrap();
        assert_eq!(validate_all.offset, 0);
        allocator.free(validate_all);
    }

    #[test]
    fn allocate_offset_allocator_reuse_trivial() {
        let mut allocator: Allocator<u32> = Allocator::new(1024 * 1024 * 256);

        // Allocator should reuse node freed by A since the allocation C fits in the same bin (using pow2 size to be sure)
        let a = allocator.allocate(1024).unwrap();
        assert_eq!(a.offset, 0);

        let b = allocator.allocate(3456).unwrap();
        assert_eq!(b.offset, 1024);

        allocator.free(a);

        let c = allocator.allocate(1024).unwrap();
        assert_eq!(c.offset, 0);

        allocator.free(c);
        allocator.free(b);

        // End: Validate that allocator has no fragmentation left. Should be 100% clean.
        let validate_all = allocator.allocate(1024 * 1024 * 256).unwrap();
        assert_eq!(validate_all.offset, 0);
        allocator.free(validate_all);
    }

    #[test]
    fn allocate_offset_allocator_reuse_complex() {
        let mut allocator: Allocator<u32> = Allocator::new(1024 * 1024 * 256);

        // Allocator should not reuse node freed by A since the allocation C doesn't fits in the same bin
        // However node D and E fit there and should reuse node from A
        let a = allocator.allocate(1024).unwrap();
        assert_eq!(a.offset, 0);

        let b = allocator.allocate(3456).unwrap();
        assert_eq!(b.offset, 1024);

        allocator.free(a);

        let c = allocator.allocate(2345).unwrap();
        assert_eq!(c.offset, 1024 + 3456);

        let d = allocator.allocate(456).unwrap();
        assert_eq!(d.offset, 0);

        let e = allocator.allocate(512).unwrap();
        assert_eq!(e.offset, 456);

        let report = allocator.storage_report();
        assert_eq!(
            report.total_free_space,
            1024 * 1024 * 256 - 3456 - 2345 - 456 - 512
        );
        assert_ne!(report.largest_free_region, report.total_free_space);

        allocator.free(c);
        allocator.free(d);
        allocator.free(b);
        allocator.free(e);

        // End: Validate that allocator has no fragmentation left. Should be 100% clean.
        let validate_all = allocator.allocate(1024 * 1024 * 256).unwrap();
        assert_eq!(validate_all.offset, 0);
        allocator.free(validate_all);
    }

    #[test]
    fn allocate_offset_allocator_zero_fragmentation() {
        let mut allocator: Allocator<u32> = Allocator::new(1024 * 1024 * 256);

        // Allocate 256x 1MB. Should fit. Then free four random slots and reallocate four slots.
        // Plus free four contiguous slots an allocate 4x larger slot. All must be zero fragmentation!
        let mut allocations: [_; 256] = std::array::from_fn(|i| {
            let allocation = allocator.allocate(1024 * 1024).unwrap();
            assert_eq!(allocation.offset, i as u32 * 1024 * 1024);
            allocation
        });

        let report = allocator.storage_report();
        assert_eq!(report.total_free_space, 0);
        assert_eq!(report.largest_free_region, 0);

        // Free four random slots
        allocator.free(allocations[243]);
        allocator.free(allocations[5]);
        allocator.free(allocations[123]);
        allocator.free(allocations[95]);

        // Free four contiguous slots (allocator must merge)
        allocator.free(allocations[151]);
        allocator.free(allocations[152]);
        allocator.free(allocations[153]);
        allocator.free(allocations[154]);

        allocations[243] = allocator.allocate(1024 * 1024).unwrap();
        allocations[5] = allocator.allocate(1024 * 1024).unwrap();
        allocations[123] = allocator.allocate(1024 * 1024).unwrap();
        allocations[95] = allocator.allocate(1024 * 1024).unwrap();
        allocations[151] = allocator.allocate(1024 * 1024 * 4).unwrap(); // 4x larger

        for (i, allocation) in allocations.iter().enumerate() {
            if !(152..155).contains(&i) {
                allocator.free(*allocation);
            }
        }

        let report2 = allocator.storage_report();
        assert_eq!(report2.total_free_space, 1024 * 1024 * 256);
        assert_eq!(report2.largest_free_region, 1024 * 1024 * 256);

        // End: Validate that allocator has no fragmentation left. Should be 100% clean.
        let validate_all = allocator.allocate(1024 * 1024 * 256).unwrap();
        assert_eq!(validate_all.offset, 0);
        allocator.free(validate_all);
    }

    #[test]
    fn min_allocator_size() {
        // Randomly generated integers on a log distribution, σ = 10.
        static TEST_OBJECT_SIZES: [u32; 42] = [
            0, 1, 2, 3, 4, 5, 8, 17, 23, 36, 51, 68, 87, 151, 165, 167, 201, 223, 306, 346, 394,
            411, 806, 969, 1404, 1798, 2236, 4281, 4745, 13989, 21095, 26594, 27146, 29679, 144685,
            153878, 495127, 727999, 1377073, 9440387, 41994490, 68520116,
        ];

        for needed_object_size in TEST_OBJECT_SIZES {
            let allocator_size = offset_allocator::min_allocator_size(needed_object_size);
            let mut allocator: Allocator<u32> = Allocator::new(allocator_size);
            assert!(allocator.allocate(needed_object_size).is_some());
        }
    }
}
