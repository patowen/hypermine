//! This module is a modified version of the https://github.com/pcwalton/offset-allocator crate,
//! updated for improved readability and optimized for Hypermine's use case.

mod bins_map;
mod node_index;
mod small_float;

use std::fmt::{Debug, Formatter};

use tracing::debug;

use crate::offset_allocator::{
    bins_map::BinsMap,
    node_index::{NodeIndex, NodeIndexOption},
    small_float::{SmallFloat, SmallFloatMap},
};

/// An allocator that manages a single contiguous chunk of space and hands out
/// portions of it as requested.
pub struct Allocator<NI: NodeIndex = u32> {
    size: u32,
    max_allocs: u32,
    /// How much available space there is across all nodes
    free_storage: u32,

    bins_map: BinsMap<NI>,

    nodes: NodeMap<NI>,
    free_nodes: FreeNodeStack<NI>,
}

struct NodeMap<NI: NodeIndex>(Vec<Node<NI>>);

impl<NI: NodeIndex> NodeMap<NI> {
    fn with_max_allocs(max_allocs: u32) -> Self {
        Self(vec![Node::default(); max_allocs as usize])
    }
}

impl<NI: NodeIndex> std::ops::Index<NI> for NodeMap<NI> {
    type Output = Node<NI>;

    fn index(&self, index: NI) -> &Self::Output {
        &self.0[index.to_usize()]
    }
}

impl<NI: NodeIndex> std::ops::IndexMut<NI> for NodeMap<NI> {
    fn index_mut(&mut self, index: NI) -> &mut Self::Output {
        &mut self.0[index.to_usize()]
    }
}

struct FreeNodeStack<NI: NodeIndex>(Vec<NI>);

impl<NI: NodeIndex> FreeNodeStack<NI> {
    fn with_max_allocs(max_allocs: u32) -> Self {
        Self((0..max_allocs).rev().map(|i| NI::from_u32(i)).collect())
    }

    #[inline]
    fn is_exhausted(&self) -> bool {
        self.0.is_empty()
    }

    #[inline]
    fn push(&mut self, node_index: NI) {
        debug!(
            "Putting node {} into freelist[{}] (free)",
            node_index,
            self.0.len()
        );
        self.0.push(node_index);
    }

    #[inline]
    fn pop_required(&mut self) -> NI {
        let node_index = self.0.pop().expect("Stack must not be exhausted");
        debug!(
            "Getting node {} from freelist[{}]",
            node_index,
            self.0.len()
        );
        node_index
    }
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
struct Node<NI = u32>
where
    NI: NodeIndex,
{
    /// The offset of the node in the address space of the heap
    data_offset: u32,
    /// The size of the node in the address space of the heap
    data_size: u32,
    /// (Patrick) Part of a linked list. Stores the previous `node_index` in the same bin
    bin_list_prev: NodeIndexOption<NI>,
    /// (Patrick) Part of a linked list. Stores the next `node_index` in the same bin
    bin_list_next: NodeIndexOption<NI>,
    /// (Patrick) Part of a linked list. Stores the previous `node_index` contiguously in memory
    neighbor_prev: NodeIndexOption<NI>,
    /// (Patrick) Part of a linked list. Stores the next `node_index` contiguously in memory
    neighbor_next: NodeIndexOption<NI>,
    /// (Patrick) Whether the node is used in an active allocation (rather than being a free slot). If `true`, this node is no longer in a bin.
    used: bool, // TODO: Merge as bit flag
}

impl<NI: NodeIndex> Default for Node<NI> {
    fn default() -> Self {
        Self {
            data_offset: Default::default(),
            data_size: Default::default(),
            bin_list_prev: Default::default(),
            bin_list_next: Default::default(),
            neighbor_prev: Default::default(),
            neighbor_next: Default::default(),
            used: Default::default(),
        }
    }
}

impl<NI> Allocator<NI>
where
    NI: NodeIndex,
{
    /// Creates a new allocator, managing a contiguous block of memory of `size`
    /// units, with a default reasonable number of maximum allocations.
    pub fn new(size: u32) -> Self {
        Allocator::with_max_allocs(size, u32::min(128 * 1024, NI::NUM_VALID))
    }

    /// Creates a new allocator, managing a contiguous block of memory of `size`
    /// units, with the given number of maximum allocations.
    ///
    /// Note that the maximum number of allocations must be less than
    /// [`NodeIndex::MAX`] minus one. If this restriction is violated, this
    /// constructor will panic.
    pub fn with_max_allocs(size: u32, max_allocs: u32) -> Self {
        assert!(max_allocs < NI::NUM_VALID);

        let mut this = Self {
            size,
            max_allocs,
            free_storage: 0,
            bins_map: BinsMap::default(),
            nodes: NodeMap::with_max_allocs(max_allocs),
            free_nodes: FreeNodeStack::with_max_allocs(max_allocs),
        };
        this.insert_node_into_bin(size, 0);
        this
    }

    /// Clears out all allocations.
    pub fn reset(&mut self) {
        *self = Self::with_max_allocs(self.size, self.max_allocs);
    }

    /// Allocates a block of `size` elements and returns its allocation.
    ///
    /// If there's not enough contiguous space for this allocation, returns
    /// None.
    pub fn allocate(&mut self, size: u32) -> Option<Allocation<NI>> {
        // Out of allocations?
        if self.free_nodes.is_exhausted() {
            // TODO: Do we want to allow an allocation that doesn't create a new node?
            return None;
        }

        // Round up to bin index to ensure that alloc >= bin
        // Gives us min bin index that fits the size
        let min_bin_index = SmallFloat::from_u32_round_up(size);
        let bin_index = self.bins_map.min_occupied_since(min_bin_index)?;

        // Pop the top node of the bin. Bin top = node.next.
        let node_index = self.bins_map[bin_index].unwrap();
        let node = &mut self.nodes[node_index];
        let node_total_size = node.data_size;
        node.data_size = size;
        node.used = true;
        self.bins_map[bin_index] = node.bin_list_next; // TODO: Doesn't this need IndexMut?
        if let Some(bin_list_next) = node.bin_list_next.to_option() {
            self.nodes[bin_list_next].bin_list_prev = NodeIndexOption::NONE;
        }
        self.free_storage -= node_total_size;
        debug!(
            "Free storage: {} (-{}) (allocate)",
            self.free_storage, node_total_size
        );

        // Bin empty?
        if self.bins_map[bin_index].is_none() {
            self.bins_map.mark_bin_empty(bin_index);
        }

        // Push back remainder N elements to a lower bin
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
    /// It may or may not panic. Note that, because this crate contains no
    /// unsafe code, the memory safe of the allocator *itself* will be
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

                // Remove node from the bin linked list and put it in the
                // freelist
                self.remove_node_from_bin(neighbor_prev);

                let prev_node = &self.nodes[neighbor_prev];
                debug_assert_eq!(prev_node.neighbor_next, NodeIndexOption::some(node_index));
                self.nodes[node_index].neighbor_prev = prev_node.neighbor_prev;
            }
        }

        if let Some(neighbor_next) = self.nodes[node_index].neighbor_next.to_option() {
            if !self.nodes[neighbor_next].used {
                // Next (contiguous) free node: Offset remains the same. Sum
                // sizes.
                let next_node = &self.nodes[neighbor_next];
                size += next_node.data_size;

                // Remove node from the bin linked list and put it in the
                // freelist
                self.remove_node_from_bin(neighbor_next);

                let next_node = &self.nodes[neighbor_next];
                debug_assert_eq!(next_node.neighbor_prev, NodeIndexOption::some(node_index));
                self.nodes[node_index].neighbor_next = next_node.neighbor_next;
            }
        }

        let Node {
            neighbor_next,
            neighbor_prev,
            ..
        } = self.nodes[node_index];

        // Insert the removed node to freelist
        self.free_nodes.push(node_index);

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

    fn insert_node_into_bin(&mut self, size: u32, data_offset: u32) -> NI {
        // Round down to bin index to ensure that bin >= alloc
        let bin_index = SmallFloat::from_u32_round_down(size);

        // Bin was empty before?
        if self.bins_map[bin_index].is_none() {
            // Set bin mask bits
            self.bins_map.mark_bin_occupied(bin_index);
        }

        // Take a freelist node and insert on top of the bin linked list (next = old top)
        let top_node_index = self.bins_map[bin_index];
        let node_index = self.free_nodes.pop_required();
        self.nodes[node_index] = Node {
            data_offset,
            data_size: size,
            bin_list_next: top_node_index,
            ..Node::default()
        };
        if let Some(top_node_index) = top_node_index.to_option() {
            self.nodes[top_node_index].bin_list_prev = NodeIndexOption::some(node_index);
        }
        self.bins_map[bin_index] = NodeIndexOption::some(node_index);

        self.free_storage += size;
        debug!(
            "Free storage: {} (+{}) (insert_node_into_bin)",
            self.free_storage, size
        );
        node_index
    }

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

                // Round down to bin index to ensure that bin >= alloc
                let bin_index = SmallFloat::from_u32_round_down(node.data_size);

                self.bins_map[bin_index] = node.bin_list_next;
                if let Some(bin_list_next) = node.bin_list_next.to_option() {
                    self.nodes[bin_list_next].bin_list_prev = NodeIndexOption::NONE;
                }

                // Bin empty?
                if self.bins_map[bin_index].is_none() {
                    self.bins_map.mark_bin_empty(bin_index);
                }
            }
        }

        // Insert the node to freelist
        self.free_nodes.push(node_index);

        self.free_storage -= node.data_size;
        debug!(
            "Free storage: {} (-{}) (remove_node_from_bin)",
            self.free_storage, node.data_size
        );
    }

    /// Returns the *used* size of an allocation.
    ///
    /// Note that this may be larger than the size requested at allocation time,
    /// due to rounding. (Patrick) No, it's never larger.
    pub fn allocation_size(&self, allocation: Allocation<NI>) -> u32 {
        self.nodes[allocation.metadata].data_size
    }

    /// Returns a structure containing the amount of free space remaining, as
    /// well as the largest amount that can be allocated at once.
    pub fn storage_report(&self) -> StorageReport {
        if self.free_nodes.is_exhausted() {
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

    /// Returns detailed information about the number of allocations in each
    /// bin.
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

impl<NI> Debug for Allocator<NI>
where
    NI: NodeIndex,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        self.storage_report().fmt(f)
    }
}

/// Returns the minimum allocator size needed to hold an object of the given
/// size.
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
