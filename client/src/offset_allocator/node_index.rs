use std::fmt::{Debug, Display};

/// The index used to identify nodes in the allocator. Determines the number of allocations
/// that the allocator supports.
///
/// By default, [`Allocator`] and related functions use `u32`, which allows for
/// `u32::MAX` allocations. You can, however, use `u16` instead, which
/// causes the allocator to use less memory but limits the number of allocations
/// within a single allocator to at most 65,535.
pub trait NodeIndex: Display + Debug + Clone + Copy + PartialEq + Eq {
    /// An invalid representation in its type, used as the `None` type of `NodeIndexOption`.
    const INVALID: Self;

    /// The number of indexes, consectuive starting from 0, that are valid representations
    const NUM_VALID: u32;

    /// Converts from a unsigned 32-bit integer to an instance of this type.
    fn from_u32(val: u32) -> Self;

    /// Converts this type to an unsigned machine word.
    fn to_usize(self) -> usize;
}

/// A type much like [`Option<NodeIndex>`] but made to use the maximum integer value as the `None`
/// value instead of requiring a separate discriminant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NodeIndexOption<NI: NodeIndex>(NI);

impl<NI: NodeIndex> NodeIndexOption<NI> {
    /// Equivalent to [`Option::None`]
    pub const NONE: Self = NodeIndexOption(NodeIndex::INVALID);

    /// Initializes what is equivalent to an [`Option::Some`] for the given node index
    pub fn some(inner: NI) -> Self {
        Self(inner)
    }

    /// Converts to the [`Option`] type for easier processing
    #[inline]
    pub fn to_option(self) -> Option<NI> {
        if self == Self::NONE {
            None
        } else {
            Some(self.0)
        }
    }

    /// Whether the option holds no value
    #[inline]
    pub fn is_none(self) -> bool {
        self == Self::NONE
    }

    /// Returns the value contained within the option. Panics if there is no such option.
    #[inline]
    pub fn unwrap(self) -> NI {
        assert!(self != Self::NONE);
        self.0
    }
}

impl<NI: NodeIndex> Default for NodeIndexOption<NI> {
    fn default() -> Self {
        Self::NONE
    }
}

impl NodeIndex for u32 {
    const INVALID: u32 = u32::MAX;

    const NUM_VALID: u32 = Self::INVALID;

    fn from_u32(val: u32) -> Self {
        assert!(val < Self::NUM_VALID);
        val
    }

    fn to_usize(self) -> usize {
        self as usize
    }
}

impl NodeIndex for u16 {
    const INVALID: u16 = u16::MAX;

    const NUM_VALID: u32 = Self::INVALID as u32;

    fn from_u32(val: u32) -> Self {
        assert!(val < Self::NUM_VALID);
        val as u16
    }

    fn to_usize(self) -> usize {
        self as usize
    }
}
