use std::{
    fmt::Debug,
    ops::{Index, IndexMut},
};

/// Like a `Vec` but with a fixed capacity
pub struct ArrayVec<T, const CAPACITY: usize> {
    contents: [T; CAPACITY],
    len: usize,
}

impl<T, const CAPACITY: usize> ArrayVec<T, CAPACITY> {
    #[inline]
    pub fn new() -> Self
    where
        T: Default + Copy,
    {
        Self::new_with_default(Default::default())
    }

    #[inline]
    pub fn new_with_default(default: T) -> Self
    where
        T: Copy,
    {
        ArrayVec {
            contents: [default; CAPACITY],
            len: 0,
        }
    }

    #[inline]
    pub fn push(&mut self, value: T) {
        self.contents[self.len] = value;
        self.len += 1;
    }

    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.contents[0..self.len].iter()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn fill(&mut self, value: T)
    where
        T: Clone,
    {
        self.contents.fill(value);
    }
}

impl<T, const CAPACITY: usize> Index<usize> for ArrayVec<T, CAPACITY> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < self.len);
        self.contents.index(index)
    }
}

impl<T, const CAPACITY: usize> IndexMut<usize> for ArrayVec<T, CAPACITY> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        assert!(index < self.len);
        self.contents.index_mut(index)
    }
}

impl<T: Debug, const CAPACITY: usize> Debug for ArrayVec<T, CAPACITY> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}
