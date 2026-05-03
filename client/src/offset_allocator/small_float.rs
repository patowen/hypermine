#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmallFloat(u32);

impl SmallFloat {
    const MANTISSA_BITS: u32 = 3;
    const EXPONENT_BITS: u32 = 5;
    const NUM_VALUES: usize = 1 << (Self::MANTISSA_BITS + Self::EXPONENT_BITS);
    const MANTISSA_VALUE: u32 = 1 << Self::MANTISSA_BITS;
    const MANTISSA_MASK: u32 = Self::MANTISSA_VALUE - 1;

    pub fn values() -> impl ExactSizeIterator<Item = Self> {
        (0..Self::NUM_VALUES).map(|i| Self(i as u32))
    }

    pub fn from_u32_round_up(value: u32) -> Self {
        let mut exp = 0;
        let mut mantissa;

        if value < Self::MANTISSA_VALUE {
            // Denorm: 0..(MANTISSA_VALUE-1)
            mantissa = value
        } else {
            // Normalized: Hidden high bit always 1. Not stored. Just like float.
            let highest_set_bit = value.ilog2();
            let mantissa_start_bit = highest_set_bit - Self::MANTISSA_BITS;
            exp = mantissa_start_bit + 1;
            mantissa = (value >> mantissa_start_bit) & Self::MANTISSA_MASK;

            let low_bits_mask = (1 << mantissa_start_bit) - 1;

            // Round up!
            if (value & low_bits_mask) != 0 {
                mantissa += 1;
            }
        }

        // + allows mantissa->exp overflow for round up
        SmallFloat((exp << Self::MANTISSA_BITS) + mantissa)
    }

    pub fn from_u32_round_down(value: u32) -> Self {
        let mut exp = 0;
        let mantissa;

        if value < Self::MANTISSA_VALUE {
            // Denorm: 0..(MANTISSA_VALUE-1)
            mantissa = value
        } else {
            // Normalized: Hidden high bit always 1. Not stored. Just like float.
            let highest_set_bit = value.ilog2();
            let mantissa_start_bit = highest_set_bit - Self::MANTISSA_BITS;
            exp = mantissa_start_bit + 1;
            mantissa = (value >> mantissa_start_bit) & Self::MANTISSA_MASK;
        }

        SmallFloat((exp << Self::MANTISSA_BITS) | mantissa)
    }

    pub fn to_u32(self) -> u32 {
        let exponent = self.0 >> Self::MANTISSA_BITS;
        let mantissa = self.0 & Self::MANTISSA_MASK;
        if exponent == 0 {
            mantissa
        } else {
            (mantissa | Self::MANTISSA_VALUE) << (exponent - 1)
        }
    }

    #[inline]
    pub fn reinterpret_as_u32(self) -> u32 {
        self.0
    }

    #[inline]
    pub fn reinterpret_u32(data: u32) -> Self {
        Self(data)
    }
}

#[derive(Debug)]
pub struct SmallFloatMap<T>([T; SmallFloat::NUM_VALUES]);

impl<T: Default + Copy> Default for SmallFloatMap<T> {
    fn default() -> Self {
        Self([T::default(); SmallFloat::NUM_VALUES])
    }
}

impl<T> std::ops::Index<SmallFloat> for SmallFloatMap<T> {
    type Output = T;

    fn index(&self, index: SmallFloat) -> &Self::Output {
        &self.0[index.0 as usize]
    }
}

impl<T> std::ops::IndexMut<SmallFloat> for SmallFloatMap<T> {
    fn index_mut(&mut self, index: SmallFloat) -> &mut Self::Output {
        &mut self.0[index.0 as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_float_uint_to_float() {
        // Denorms, exp=1 and exp=2 + mantissa = 0 are all precise.
        // NOTE: Assuming 8 value (3 bit) mantissa.
        // If this test fails, please change this assumption!
        let precise_number_count = 17;
        for i in 0..precise_number_count {
            let round_up = SmallFloat::from_u32_round_up(i);
            let round_down = SmallFloat::from_u32_round_down(i);
            assert_eq!(SmallFloat::reinterpret_u32(i), round_up);
            assert_eq!(SmallFloat::reinterpret_u32(i), round_down);
        }

        // Test some random picked numbers
        struct NumberFloatUpDown {
            number: u32,
            up: SmallFloat,
            down: SmallFloat,
        }

        let test_data = [
            NumberFloatUpDown {
                number: 17,
                up: SmallFloat::reinterpret_u32(17),
                down: SmallFloat::reinterpret_u32(16),
            },
            NumberFloatUpDown {
                number: 118,
                up: SmallFloat::reinterpret_u32(39),
                down: SmallFloat::reinterpret_u32(38),
            },
            NumberFloatUpDown {
                number: 1024,
                up: SmallFloat::reinterpret_u32(64),
                down: SmallFloat::reinterpret_u32(64),
            },
            NumberFloatUpDown {
                number: 65536,
                up: SmallFloat::reinterpret_u32(112),
                down: SmallFloat::reinterpret_u32(112),
            },
            NumberFloatUpDown {
                number: 529445,
                up: SmallFloat::reinterpret_u32(137),
                down: SmallFloat::reinterpret_u32(136),
            },
            NumberFloatUpDown {
                number: 1048575,
                up: SmallFloat::reinterpret_u32(144),
                down: SmallFloat::reinterpret_u32(143),
            },
        ];

        for v in test_data {
            let round_up = SmallFloat::from_u32_round_up(v.number);
            let round_down = SmallFloat::from_u32_round_down(v.number);
            assert_eq!(round_up, v.up);
            assert_eq!(round_down, v.down);
        }
    }

    #[test]
    fn small_float_float_to_uint() {
        // Denorms, exp=1 and exp=2 + mantissa = 0 are all precise.
        // NOTE: Assuming 8 value (3 bit) mantissa.
        // If this test fails, please change this assumption!
        let precise_number_count = 17;
        for i in 0..precise_number_count {
            let v = SmallFloat::reinterpret_u32(i).to_u32();
            assert_eq!(i, v);
        }

        // Test that float->uint->float conversion is precise for all numbers
        // NOTE: Test values < 240. 240->4G = overflows 32 bit integer
        for i in (0..240).map(|i| SmallFloat::reinterpret_u32(i)) {
            let v = i.to_u32();
            let round_up = SmallFloat::from_u32_round_up(v);
            let round_down = SmallFloat::from_u32_round_down(v);
            assert_eq!(i, round_up);
            assert_eq!(i, round_down);
        }
    }
}
