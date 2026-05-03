use crate::offset_allocator::{
    node_index::{NodeIndex, NodeIndexOption},
    small_float::{SmallFloat, SmallFloatMap},
};

const NUM_TOP_BINS: usize = 32;
const TOP_BINS_INDEX_SHIFT: u32 = 3;
const LEAF_BINS_INDEX_MASK: u32 = 7;

/// A map from each bin to the node at the head of the linked list for that bin. The name of this struct is `BinsMap` instead of `BinMap` to avoid confusion with binary maps.
pub struct BinsMap<NI: NodeIndex> {
    /// (Patrick) A bit-vector showing which `top_bin_index`es are "used"
    occupied_bins_top: u32,
    /// (Patrick) An array of 32 bit-vectors showing which `leaf_bin_index`es are "used" for the given top bin, usually indexed by `top_bin_index`
    occupied_bins: [u8; NUM_TOP_BINS],
    /// (Patrick) An array of 256 `node_index`es (each being a head of a doubly-linked list of nodes in that bin), usually indexed by `bin_index` (a combo of `top_bin_index` and `leaf_bin_index`).
    head_nodes: SmallFloatMap<NodeIndexOption<NI>>,
}

impl<NI: NodeIndex> Default for BinsMap<NI> {
    fn default() -> Self {
        Self {
            occupied_bins_top: 0,
            occupied_bins: [0; NUM_TOP_BINS],
            head_nodes: SmallFloatMap::default(),
        }
    }
}

impl<NI: NodeIndex> BinsMap<NI> {
    pub fn min_occupied_since(&self, min: SmallFloat) -> Option<SmallFloat> {
        /// Find the lowest bit that is set to 1, as long as it's at least start_bit_index. Return `None` if there is no such bit.
        fn find_lowest_bit_set_after(bit_mask: u32, start_bit_index: u32) -> Option<u32> {
            let mask_before_start_index = (1 << start_bit_index) - 1;
            let mask_after_start_index = !mask_before_start_index;
            let bits_after = bit_mask & mask_after_start_index;
            if bits_after == 0 {
                None
            } else {
                Some(bits_after.trailing_zeros())
            }
        }

        let min_top_bin_index = min.reinterpret_as_u32() >> TOP_BINS_INDEX_SHIFT;
        let min_leaf_bin_index = min.reinterpret_as_u32() & LEAF_BINS_INDEX_MASK;

        let mut top_bin_index = min_top_bin_index;
        let mut leaf_bin_index = None;

        // If top bin exists, scan its leaf bin. This can fail (NO_SPACE).
        if (self.occupied_bins_top & (1 << top_bin_index)) != 0 {
            leaf_bin_index = find_lowest_bit_set_after(
                self.occupied_bins[top_bin_index as usize] as _,
                min_leaf_bin_index,
            );
        }

        // If we didn't find space in top bin, we search top bin from +1
        let leaf_bin_index = match leaf_bin_index {
            Some(leaf_bin_index) => leaf_bin_index,
            None => {
                top_bin_index =
                    find_lowest_bit_set_after(self.occupied_bins_top, min_top_bin_index + 1)?;

                // All leaf bins here fit the alloc, since the top bin was
                // rounded up. Start leaf search from bit 0.
                //
                // NOTE: This search can't fail since at least one leaf bit was
                // set because the top bit was set.
                self.occupied_bins[top_bin_index as usize].trailing_zeros()
            }
        };

        Some(SmallFloat::reinterpret_u32(
            (top_bin_index << TOP_BINS_INDEX_SHIFT) | leaf_bin_index,
        ))
    }

    pub fn max_occupied(&self) -> Option<SmallFloat> {
        if self.occupied_bins_top == 0 {
            return None;
        }
        let top_bin_index = self.occupied_bins_top.ilog2();
        let leaf_bin_index = (self.occupied_bins[top_bin_index as usize] as u32).ilog2();
        Some(SmallFloat::reinterpret_u32(
            (top_bin_index << TOP_BINS_INDEX_SHIFT) | leaf_bin_index,
        ))
    }

    pub fn mark_bin_empty(&mut self, bin_index: SmallFloat) {
        let top_bin_index = bin_index.reinterpret_as_u32() >> TOP_BINS_INDEX_SHIFT;
        let leaf_bin_index = bin_index.reinterpret_as_u32() & LEAF_BINS_INDEX_MASK;

        // Remove a leaf bin mask bit
        self.occupied_bins[top_bin_index as usize] &= !(1 << u32::from(leaf_bin_index));

        // All leaf bins empty?
        if self.occupied_bins[top_bin_index as usize] == 0 {
            // Remove a top bin mask bit
            self.occupied_bins_top &= !(1 << top_bin_index);
        }
    }

    pub fn mark_bin_occupied(&mut self, bin_index: SmallFloat) {
        let top_bin_index = bin_index.reinterpret_as_u32() >> TOP_BINS_INDEX_SHIFT;
        let leaf_bin_index = bin_index.reinterpret_as_u32() & LEAF_BINS_INDEX_MASK;

        // Set bin mask bits
        self.occupied_bins[top_bin_index as usize] |= 1 << leaf_bin_index;
        self.occupied_bins_top |= 1 << top_bin_index;
    }
}

impl<NI: NodeIndex> std::ops::Index<SmallFloat> for BinsMap<NI> {
    type Output = NodeIndexOption<NI>;

    fn index(&self, index: SmallFloat) -> &Self::Output {
        &self.head_nodes[index]
    }
}

impl<NI: NodeIndex> std::ops::IndexMut<SmallFloat> for BinsMap<NI> {
    fn index_mut(&mut self, index: SmallFloat) -> &mut Self::Output {
        &mut self.head_nodes[index]
    }
}
