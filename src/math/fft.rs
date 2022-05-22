//! The Fast Fourier Transform (FFT) and Number Theoretic Transform (NTT)
use super::num::{CommonField, Complex, PI};
use std::ops::{Add, Div, Mul, Neg, Sub};

// We can delete this struct once f64::reverse_bits() stabilizes.
struct BitRevIterator {
    a: usize,
    n: usize,
}
impl BitRevIterator {
    fn new(n: usize) -> Self {
        assert!(n.is_power_of_two());
        Self { a: 2 * n - 1, n }
    }
}
impl Iterator for BitRevIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.a == 2 * self.n - 2 {
            return None;
        }
        let mut mask = self.n;
        while self.a & mask > 0 {
            self.a ^= mask;
            mask /= 2;
        }
        self.a |= mask;
        Some(self.a / 2)
    }
}

#[allow(clippy::upper_case_acronyms)]
pub trait FFT: Sized + Copy {
    type F: Sized
        + Copy
        + From<Self>
        + Neg
        + Add<Output = Self::F>
        + Div<Output = Self::F>
        + Mul<Output = Self::F>
        + Sub<Output = Self::F>;

    const ZERO: Self;

    /// A primitive nth root of one raised to the powers 0, 1, 2, ..., n/2 - 1
    fn get_roots(n: usize, inverse: bool) -> Vec<Self::F>;
    /// 1 for forward transform, 1/n for inverse transform
    fn get_factor(n: usize, inverse: bool) -> Self::F;
    /// The inverse of Self::F::from()
    fn extract(f: Self::F) -> Self;
}

impl FFT for f64 {
    type F = Complex;

    const ZERO: f64 = 0.0;

    fn get_roots(n: usize, inverse: bool) -> Vec<Self::F> {
        let step = if inverse { -2.0 } else { 2.0 } * PI / n as f64;
        (0..n / 2)
            .map(|i| Complex::from_polar(1.0, step * i as f64))
            .collect()
    }

    fn get_factor(n: usize, inverse: bool) -> Self::F {
        Self::F::from(if inverse { (n as f64).recip() } else { 1.0 })
    }

    fn extract(f: Self::F) -> f64 {
        f.real
    }
}

// NTT notes: see problem 30-6 in CLRS for details, keeping in mind that
//      2187 and  410692747 are inverses and 2^26th roots of 1 mod (7<<26)+1
//  15311432 and  469870224 are inverses and 2^23rd roots of 1 mod (119<<23)+1
// 440564289 and 1713844692 are inverses and 2^27th roots of 1 mod (15<<27)+1
//       125 and 2267742733 are inverses and 2^30th roots of 1 mod (3<<30)+1
impl FFT for i64 {
    type F = CommonField;

    const ZERO: Self = 0;

    fn get_roots(n: usize, inverse: bool) -> Vec<Self::F> {
        assert!(n <= 1 << 23);
        let mut prim_root = Self::F::from(15_311_432);
        if inverse {
            prim_root = prim_root.recip();
        }
        for _ in (0..).take_while(|&i| n < 1 << (23 - i)) {
            prim_root = prim_root * prim_root;
        }

        let mut roots = Vec::with_capacity(n / 2);
        let mut root = Self::F::from(1);
        for _ in 0..roots.capacity() {
            roots.push(root);
            root = root * prim_root;
        }
        roots
    }

    fn get_factor(n: usize, inverse: bool) -> Self::F {
        Self::F::from(if inverse { n as Self } else { 1 }).recip()
    }

    fn extract(f: Self::F) -> Self {
        f.val
    }
}

/// Computes the discrete fourier transform of v, whose length is a power of 2.
/// Forward transform: polynomial coefficients -> evaluate at roots of unity
/// Inverse transform: values at roots of unity -> interpolated coefficients
pub fn fft<T: FFT>(v: &[T::F], inverse: bool) -> Vec<T::F> {
    let n = v.len();
    assert!(n.is_power_of_two());

    let factor = T::get_factor(n, inverse);
    let roots_of_unity = T::get_roots(n, inverse);
    let mut dft = BitRevIterator::new(n)
        .map(|i| v[i] * factor)
        .collect::<Vec<_>>();

    for m in (0..).map(|s| 1 << s).take_while(|&m| m < n) {
        for k in (0..n).step_by(2 * m) {
            for j in 0..m {
                let u = dft[k + j];
                let t = dft[k + j + m] * roots_of_unity[n / 2 / m * j];
                dft[k + j] = u + t;
                dft[k + j + m] = u - t;
            }
        }
    }
    dft
}

/// From a slice of reals (f64 or i64), computes DFT of size at least desired_len
pub fn dft_from_reals<T: FFT>(v: &[T], desired_len: usize) -> Vec<T::F> {
    assert!(v.len() <= desired_len);

    let complex_v = v
        .iter()
        .cloned()
        .chain(std::iter::repeat(T::ZERO))
        .take(desired_len.next_power_of_two())
        .map(T::F::from)
        .collect::<Vec<_>>();
    fft::<T>(&complex_v, false)
}

/// The inverse of dft_from_reals()
pub fn idft_to_reals<T: FFT>(dft_v: &[T::F], desired_len: usize) -> Vec<T> {
    assert!(dft_v.len() >= desired_len);

    let complex_v = fft::<T>(dft_v, true);
    complex_v
        .into_iter()
        .take(desired_len)
        .map(T::extract)
        .collect()
}

/// Given two polynomials (vectors) sum_i a[i] x^i and sum_i b[i] x^i,
/// computes their product (convolution) c[k] = sum_(i+j=k) a[i]*b[j].
/// Uses complex FFT if inputs are f64, or modular NTT if inputs are i64.
pub fn convolution<T: FFT>(a: &[T], b: &[T]) -> Vec<T> {
    let len_c = a.len() + b.len() - 1;
    let dft_a = dft_from_reals(a, len_c).into_iter();
    let dft_b = dft_from_reals(b, len_c).into_iter();
    let dft_c = dft_a.zip(dft_b).map(|(a, b)| a * b).collect::<Vec<_>>();
    idft_to_reals(&dft_c, len_c)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_complex_dft() {
        let v = vec![7.0, 1.0, 1.0];
        let dft_v = dft_from_reals(&v, v.len());
        let new_v: Vec<f64> = idft_to_reals(&dft_v, v.len());

        let six = Complex::from(6.0);
        let seven = Complex::from(7.0);
        let nine = Complex::from(9.0);
        let i = Complex::new(0.0, 1.0);

        assert_eq!(dft_v, vec![nine, six + i, seven, six - i]);
        assert_eq!(new_v, v);
    }

    #[test]
    fn test_modular_dft() {
        let v = vec![7, 1, 1];
        let dft_v = dft_from_reals(&v, v.len());
        let new_v: Vec<i64> = idft_to_reals(&dft_v, v.len());

        let seven = CommonField::from(7);
        let one = CommonField::from(1);
        let prim = CommonField::from(15_311_432).pow(1 << 21);
        let prim2 = prim * prim;

        let eval0 = seven + one + one;
        let eval1 = seven + prim + prim2;
        let eval2 = seven + prim2 + one;
        let eval3 = seven + prim.recip() + prim2;

        assert_eq!(dft_v, vec![eval0, eval1, eval2, eval3]);
        assert_eq!(new_v, v);
    }

    #[test]
    fn test_complex_convolution() {
        let x = vec![7.0, 1.0, 1.0];
        let y = vec![2.0, 4.0];
        let z = convolution(&x, &y);
        let m = convolution(&vec![999.0], &vec![1e6]);

        assert_eq!(z, vec![14.0, 30.0, 6.0, 4.0]);
        assert_eq!(m, vec![999e6]);
    }

    #[test]
    fn test_modular_convolution() {
        let x = vec![7, 1, 1];
        let y = vec![2, 4];
        let z = convolution(&x, &y);
        let m = convolution(&vec![999], &vec![1_000_000]);

        assert_eq!(z, vec![14, 30, 6, 4]);
        assert_eq!(m, vec![999_000_000 - super::super::num::COMMON_PRIME]);
    }
}
