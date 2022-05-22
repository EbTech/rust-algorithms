//! Rational and Complex numbers, safe modular arithmetic, and linear algebra,
//! implemented minimally for contest use.
//! If you need more features, you might be interested in crates.io/crates/num
pub use std::f64::consts::PI;
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

/// Fast iterative version of Euclid's GCD algorithm
pub fn fast_gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        a %= b;
        std::mem::swap(&mut a, &mut b);
    }
    a.abs()
}

/// Represents a fraction reduced to lowest terms
#[derive(Clone, Copy, Eq, PartialEq, Debug, Hash)]
pub struct Rational {
    pub num: i64,
    pub den: i64,
}
impl Rational {
    pub fn new(num: i64, den: i64) -> Self {
        let g = fast_gcd(num, den) * den.signum();
        Self {
            num: num / g,
            den: den / g,
        }
    }
    pub fn abs(self) -> Self {
        Self {
            num: self.num.abs(),
            den: self.den,
        }
    }
    pub fn recip(self) -> Self {
        let g = self.num.signum();
        Self {
            num: self.den / g,
            den: self.num / g,
        }
    }
}
impl From<i64> for Rational {
    fn from(num: i64) -> Self {
        Self { num, den: 1 }
    }
}
impl Neg for Rational {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            num: -self.num,
            den: self.den,
        }
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Add for Rational {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::new(
            self.num * other.den + self.den * other.num,
            self.den * other.den,
        )
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Sub for Rational {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self::new(
            self.num * other.den - self.den * other.num,
            self.den * other.den,
        )
    }
}
impl Mul for Rational {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self::new(self.num * other.num, self.den * other.den)
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for Rational {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        self * other.recip()
    }
}
impl Ord for Rational {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.num * other.den).cmp(&(self.den * other.num))
    }
}
impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Represents a complex number using floating-point arithmetic
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Complex {
    pub real: f64,
    pub imag: f64,
}
impl Complex {
    pub fn new(real: f64, imag: f64) -> Self {
        Self { real, imag }
    }
    pub fn from_polar(r: f64, th: f64) -> Self {
        Self::new(r * th.cos(), r * th.sin())
    }
    pub fn abs_square(self) -> f64 {
        self.real * self.real + self.imag * self.imag
    }
    pub fn argument(self) -> f64 {
        self.imag.atan2(self.real)
    }
    pub fn conjugate(self) -> Self {
        Self::new(self.real, -self.imag)
    }
    pub fn recip(self) -> Self {
        let denom = self.abs_square();
        Self::new(self.real / denom, -self.imag / denom)
    }
}
impl From<f64> for Complex {
    fn from(real: f64) -> Self {
        Self::new(real, 0.0)
    }
}
impl Neg for Complex {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.real, -self.imag)
    }
}
impl Add for Complex {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::new(self.real + other.real, self.imag + other.imag)
    }
}
impl Sub for Complex {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self::new(self.real - other.real, self.imag - other.imag)
    }
}
impl Mul for Complex {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let real = self.real * other.real - self.imag * other.imag;
        let imag = self.imag * other.real + self.real * other.imag;
        Self::new(real, imag)
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for Complex {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        self * other.recip()
    }
}

/// Represents an element of the finite (Galois) field of prime order M, where
/// 1 <= M < 2^31.5. If M is not prime, ring operations are still valid
/// but recip() and division are not. Note that the latter operations are also
/// the slowest, so precompute any inverses that you intend to use frequently.
#[derive(Clone, Copy, Eq, PartialEq, Debug, Hash)]
pub struct Modulo<const M: i64> {
    pub val: i64,
}
impl<const M: i64> Modulo<M> {
    /// Computes self^n in O(log n) time
    pub fn pow(mut self, mut n: u64) -> Self {
        let mut result = Self::from_small(1);
        while n > 0 {
            if n % 2 == 1 {
                result = result * self;
            }
            self = self * self;
            n /= 2;
        }
        result
    }
    /// Computes inverses of 1 to n in O(n) time
    pub fn vec_of_recips(n: i64) -> Vec<Self> {
        let mut recips = vec![Self::from(0), Self::from(1)];
        for i in 2..=n {
            let (md, dv) = (M % i, M / i);
            recips.push(recips[md as usize] * Self::from_small(-dv));
        }
        recips
    }
    /// Computes self^-1 in O(log M) time
    pub fn recip(self) -> Self {
        self.pow(M as u64 - 2)
    }
    /// Avoids the % operation but requires -M <= x < M
    fn from_small(s: i64) -> Self {
        let val = if s < 0 { s + M } else { s };
        Self { val }
    }
}
impl<const M: i64> From<i64> for Modulo<M> {
    fn from(val: i64) -> Self {
        // Self { val: val.rem_euclid(M) }
        Self::from_small(val % M)
    }
}
impl<const M: i64> Neg for Modulo<M> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::from_small(-self.val)
    }
}
impl<const M: i64> Add for Modulo<M> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::from_small(self.val + other.val - M)
    }
}
impl<const M: i64> Sub for Modulo<M> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self::from_small(self.val - other.val)
    }
}
impl<const M: i64> Mul for Modulo<M> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self::from(self.val * other.val)
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl<const M: i64> Div for Modulo<M> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        self * other.recip()
    }
}

/// Prime modulus that's commonly used in programming competitions
pub const COMMON_PRIME: i64 = 998_244_353; // 2^23 * 7 * 17 + 1;
pub type CommonField = Modulo<COMMON_PRIME>;

#[derive(Clone, PartialEq, Debug)]
pub struct Matrix {
    cols: usize,
    inner: Box<[f64]>,
}
impl Matrix {
    pub fn zero(rows: usize, cols: usize) -> Self {
        let inner = vec![0.0; rows * cols].into_boxed_slice();
        Self { cols, inner }
    }
    pub fn one(cols: usize) -> Self {
        let mut matrix = Self::zero(cols, cols);
        for i in 0..cols {
            matrix[i][i] = 1.0;
        }
        matrix
    }
    pub fn vector(vec: &[f64], as_row: bool) -> Self {
        let cols = if as_row { vec.len() } else { 1 };
        let inner = vec.to_vec().into_boxed_slice();
        Self { cols, inner }
    }
    pub fn pow(&self, mut n: u64) -> Self {
        let mut base = self.clone();
        let mut result = Self::one(self.cols);
        while n > 0 {
            if n % 2 == 1 {
                result = &result * &base;
            }
            base = &base * &base;
            n /= 2;
        }
        result
    }
    pub fn rows(&self) -> usize {
        self.inner.len() / self.cols
    }
    pub fn transpose(&self) -> Self {
        let mut matrix = Matrix::zero(self.cols, self.rows());
        for i in 0..self.rows() {
            for j in 0..self.cols {
                matrix[j][i] = self[i][j];
            }
        }
        matrix
    }
    pub fn recip(&self) -> Self {
        unimplemented!();
    }
}
impl Index<usize> for Matrix {
    type Output = [f64];
    fn index(&self, row: usize) -> &Self::Output {
        let start = self.cols * row;
        &self.inner[start..start + self.cols]
    }
}
impl IndexMut<usize> for Matrix {
    fn index_mut(&mut self, row: usize) -> &mut Self::Output {
        let start = self.cols * row;
        &mut self.inner[start..start + self.cols]
    }
}
impl Neg for &Matrix {
    type Output = Matrix;
    fn neg(self) -> Matrix {
        let inner = self.inner.iter().map(|&v| -v).collect();
        Matrix {
            cols: self.cols,
            inner,
        }
    }
}
impl Add for &Matrix {
    type Output = Matrix;
    fn add(self, other: Self) -> Matrix {
        let self_iter = self.inner.iter();
        let inner = self_iter
            .zip(other.inner.iter())
            .map(|(&u, &v)| u + v)
            .collect();
        Matrix {
            cols: self.cols,
            inner,
        }
    }
}
impl Sub for &Matrix {
    type Output = Matrix;
    fn sub(self, other: Self) -> Matrix {
        let self_iter = self.inner.iter();
        let inner = self_iter
            .zip(other.inner.iter())
            .map(|(&u, &v)| u - v)
            .collect();
        Matrix {
            cols: self.cols,
            inner,
        }
    }
}
impl Mul<f64> for &Matrix {
    type Output = Matrix;
    fn mul(self, scalar: f64) -> Matrix {
        let inner = self.inner.iter().map(|&v| v * scalar).collect();
        Matrix {
            cols: self.cols,
            inner,
        }
    }
}
impl Mul for &Matrix {
    type Output = Matrix;
    fn mul(self, other: Self) -> Matrix {
        assert_eq!(self.cols, other.rows());
        let mut matrix = Matrix::zero(self.rows(), other.cols);
        for i in 0..self.rows() {
            for k in 0..self.cols {
                for j in 0..other.cols {
                    matrix[i][j] += self[i][k] * other[k][j];
                }
            }
        }
        matrix
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_rational() {
        let three = Rational::from(3);
        let six = Rational::from(6);
        let three_and_half = three + three / six;

        assert_eq!(three_and_half.num, 7);
        assert_eq!(three_and_half.den, 2);
        assert_eq!(three_and_half, Rational::new(-35, -10));
        assert!(three_and_half > Rational::from(3));
        assert!(three_and_half < Rational::from(4));

        let minus_three_and_half = six - three_and_half + three / (-three / six);
        let zero = three_and_half + minus_three_and_half;

        assert_eq!(minus_three_and_half.num, -7);
        assert_eq!(minus_three_and_half.den, 2);
        assert_eq!(three_and_half, -minus_three_and_half);
        assert_eq!(zero.num, 0);
        assert_eq!(zero.den, 1);
    }

    #[test]
    fn test_complex() {
        let four = Complex::new(4.0, 0.0);
        let two_i = Complex::new(0.0, 2.0);

        assert_eq!(four / two_i, -two_i);
        assert_eq!(two_i * -two_i, four);
        assert_eq!(two_i - two_i, Complex::from(0.0));
        assert_eq!(four.abs_square(), 16.0);
        assert_eq!(two_i.abs_square(), 4.0);
        assert_eq!((-four).argument(), -PI);
        assert_eq!((-two_i).argument(), -PI / 2.0);
        assert_eq!(four.argument(), 0.0);
        assert_eq!(two_i.argument(), PI / 2.0);
    }

    #[test]
    fn test_field() {
        let base = CommonField::from(1234);
        let zero = base - base;
        let one = base.recip() * base;
        let two = CommonField::from(2 - 5 * COMMON_PRIME);

        assert_eq!(zero.val, 0);
        assert_eq!(one.val, 1);
        assert_eq!(one + one, two);
        assert_eq!(one / base * (base * base) - base / one, zero);
    }

    #[test]
    fn test_vec_of_recips() {
        let recips = CommonField::vec_of_recips(20);

        assert_eq!(recips.len(), 21);
        for i in 1..recips.len() {
            assert_eq!(recips[i], CommonField::from(i as i64).recip());
        }
    }

    #[test]
    fn test_linalg() {
        let zero = Matrix::zero(2, 2);
        let one = Matrix::one(2);
        let rotate_90 = Matrix {
            cols: 2,
            inner: Box::new([0.0, -1.0, 1.0, 0.0]),
        };
        let x_vec = Matrix::vector(&[1.0, 0.0], false);
        let y_vec = Matrix::vector(&[0.0, 1.0], false);
        let x_dot_x = &x_vec.transpose() * &x_vec;
        let x_dot_y = &x_vec.transpose() * &y_vec;

        assert_eq!(x_dot_x, Matrix::one(1));
        assert_eq!(x_dot_x[0][0], 1.0);
        assert_eq!(x_dot_y, Matrix::zero(1, 1));
        assert_eq!(x_dot_y[0][0], 0.0);
        assert_eq!(&one - &one, zero);
        assert_eq!(&one * 0.0, zero);
        assert_eq!(&rotate_90 * &rotate_90, -&one);
        assert_eq!(&rotate_90 * &x_vec, y_vec);
        assert_eq!(&rotate_90 * &y_vec, -&x_vec);
        assert_eq!(&rotate_90 * &(&x_vec + &y_vec), &y_vec - &x_vec);
    }
}
