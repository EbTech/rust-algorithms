//! The Fast Fourier Transform (FFT)
use super::num::{Complex, PI};

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

// Integer FFT notes: see problem 30-6 in CLRS for details, noting that
// 440564289 and 1713844692 are inverses and 2^27th roots of 1 mod p=(15<<27)+1
//       125 and 2267742733 are inverses and 2^30th roots of 1 mod p=(3<<30)+1

/// Computes the discrete fourier transform of v, whose length is a power of 2.
/// Forward transform: polynomial coefficients -> evaluate at roots of unity
/// Inverse transform: values at roots of unity -> interpolated coefficients
pub fn fft(v: &[Complex], inverse: bool) -> Vec<Complex> {
    let n = v.len();
    assert!(n.is_power_of_two());

    let step = if inverse { -2.0 } else { 2.0 } * PI / n as f64;
    let factor = Complex::from(if inverse { n as f64 } else { 1.0 });
    let roots_of_unity = (0..n / 2)
        .map(|i| Complex::from_polar(1.0, step * i as f64))
        .collect::<Vec<_>>();
    let mut dft = BitRevIterator::new(n)
        .map(|i| v[i] / factor)
        .collect::<Vec<_>>();

    for m in (0..).map(|s| 1 << s).take_while(|&m| m < n) {
        for j in 0..m {
            for k in (0..n).step_by(2 * m) {
                let u = dft[k + j];
                let t = dft[k + j + m] * roots_of_unity[n / 2 / m * j];
                dft[k + j] = u + t;
                dft[k + j + m] = u - t;
            }
        }
    }
    dft
}

/// From a real vector, computes a DFT of size at least desired_len
pub fn dft_from_reals(v: &[f64], desired_len: usize) -> Vec<Complex> {
    assert!(v.len() <= desired_len);
    let complex_v = v
        .iter()
        .cloned()
        .chain(std::iter::repeat(0.0))
        .take(desired_len.next_power_of_two())
        .map(Complex::from)
        .collect::<Vec<_>>();
    fft(&complex_v, false)
}

/// The inverse of dft_from_reals()
pub fn idft_to_reals(fft_v: &[Complex], desired_len: usize) -> Vec<f64> {
    assert!(fft_v.len() >= desired_len);
    let complex_v = fft(fft_v, true);
    complex_v
        .into_iter()
        .take(desired_len)
        .map(|c| c.real) // to get integers: c.real.round() as i64
        .collect()
}

/// Given two polynomials (vectors) sum_i a[i]x^i and sum_i b[i]x^i,
/// computes their product (convolution) c[k] = sum_(i+j=k) a[i]*b[j]
pub fn convolution(a: &[f64], b: &[f64]) -> Vec<f64> {
    let len_c = a.len() + b.len() - 1;
    let dft_a = dft_from_reals(a, len_c);
    let dft_b = dft_from_reals(b, len_c);
    let dft_c = dft_a
        .into_iter()
        .zip(dft_b.into_iter())
        .map(|(a, b)| a * b)
        .collect::<Vec<_>>();
    idft_to_reals(&dft_c, len_c)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_dft() {
        let v = vec![7.0, 1.0, 1.0];
        let dft_v = dft_from_reals(&v, v.len());
        let new_v = idft_to_reals(&dft_v, v.len());

        let six = Complex::from(6.0);
        let seven = Complex::from(7.0);
        let nine = Complex::from(9.0);
        let i = Complex::new(0.0, 1.0);

        assert_eq!(dft_v, vec![nine, six + i, seven, six - i]);
        assert_eq!(new_v, v);
    }

    #[test]
    fn test_convolution() {
        let x = vec![2.0, 3.0, 2.0];
        let y = vec![7.0, 2.0];
        let z = convolution(&x, &y);

        assert_eq!(z, vec![14.0, 25.0, 20.0, 4.0]);
    }
}
