//! Pseudorandom number generators (PRNGs).

/// A simple and efficient random number generator.
pub type SmallRng = Xoshiro256PlusPlus;

/// A xoshiro256++ random number generator.
///
/// This is a simplified version of the `SmallRng` implementation from the
/// excellent `rand` crate, keeping only essential features.
///
/// The xoshiro256++ algorithm is not suitable for cryptographic purposes, but
/// is very fast and has excellent statistical properties.
///
/// * Source: [Docs.rs](https://docs.rs/rand/0.8.4/src/rand/rngs/xoshiro256plusplus.rs.html)
/// * Theory: [Xorshift - Wikipedia](https://en.wikipedia.org/wiki/Xorshift)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Xoshiro256PlusPlus {
    s: [u64; 4],
}

impl Xoshiro256PlusPlus {
    /// Construct a new RNG from a 64-bit seed.
    pub fn new(mut state: u64) -> Self {
        const PHI: u64 = 0x9e3779b97f4a7c15;
        let mut seed = <[u64; 4]>::default();
        for chunk in &mut seed {
            state = state.wrapping_add(PHI);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            z = z ^ (z >> 31);
            *chunk = z;
        }
        Self { s: seed }
    }

    /// Generate a random `u32`.
    #[inline]
    pub fn next_u32(&mut self) -> u32 {
        (self.next_u64() >> 32) as u32
    }

    /// Generate a random `u64`.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result_plusplus = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);

        let t = self.s[1] << 17;

        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];

        self.s[2] ^= t;

        self.s[3] = self.s[3].rotate_left(45);

        result_plusplus
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xoshiro256plusplus() {
        let mut rng = Xoshiro256PlusPlus::new(42);
        assert_eq!(rng.next_u64(), 15021278609987233951);
        assert_eq!(rng.next_u64(), 5881210131331364753);
        assert_eq!(rng.next_u64(), 18149643915985481100);
        assert_eq!(rng.next_u64(), 12933668939759105464);
        assert_eq!(rng.next_u64(), 14637574242682825331);
        assert_eq!(rng.next_u64(), 10848501901068131965);
        assert_eq!(rng.next_u64(), 2312344417745909078);
        assert_eq!(rng.next_u64(), 11162538943635311430);
    }

    #[test]
    fn reference() {
        let mut rng = Xoshiro256PlusPlus { s: [1, 2, 3, 4] };
        // These values were produced with the reference implementation:
        // http://xoshiro.di.unimi.it/xoshiro256plusplus.c
        let expected = [
            41943041,
            58720359,
            3588806011781223,
            3591011842654386,
            9228616714210784205,
            9973669472204895162,
            14011001112246962877,
            12406186145184390807,
            15849039046786891736,
            10450023813501588000,
        ];
        for &e in &expected {
            assert_eq!(rng.next_u64(), e);
        }
    }
}
