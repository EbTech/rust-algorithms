//! Safe modular arithmetic as well as Rational and Complex numbers,
//! implemented minimally for contest use.
//! If you need more features, you might be interested in crates.io/crates/num
pub use std::f64::consts::PI;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Fast iterative version of Euclid's GCD algorithm
pub fn fast_gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        a %= b;
        std::mem::swap(&mut a, &mut b);
    }
    a.abs()
}

/// Represents a fraction reduced to lowest terms
#[derive(Clone, Copy, Eq, PartialEq, Debug)]
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
    fn recip(self) -> Self {
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

/// Represents an element of the finite (Galois) field of prime order, given by
/// MOD. Until Rust gets const generics, MOD must be hardcoded, but any prime
/// in [1, 2^32] will work. If MOD is not prime, ring operations are still valid
/// but recip() and division are not. Note that the latter operations are also
/// the slowest, so precompute any inverses that you intend to use frequently.
#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub struct Field {
    pub val: u64,
}
impl Field {
    pub const MOD: u64 = 998_244_353; // 2^23 * 7 * 17 + 1

    pub fn pow(mut self, mut exp: u64) -> Self {
        let mut result = Self::from_small(1);
        while exp > 0 {
            if exp % 2 == 1 {
                result = result * self;
            }
            self = self * self;
            exp /= 2;
        }
        result
    }
    pub fn recip(self) -> Self {
        self.pow(Self::MOD - 2)
    }
    fn from_small(s: u64) -> Self {
        let val = if s < Self::MOD { s } else { s - Self::MOD };
        Self { val }
    }
}
impl From<u64> for Field {
    fn from(val: u64) -> Self {
        Self {
            val: val % Self::MOD,
        }
    }
}
impl Neg for Field {
    type Output = Self;
    fn neg(self) -> Self {
        Self::from_small(Self::MOD - self.val)
    }
}
impl Add for Field {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::from_small(self.val + other.val)
    }
}
impl Sub for Field {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self::from_small(self.val + Self::MOD - other.val)
    }
}
impl Mul for Field {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self::from(self.val * other.val)
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for Field {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        self * other.recip()
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
        let base = Field::from(1234);
        let zero = base - base;
        let one = base.recip() * base;

        assert_eq!(zero.val, 0);
        assert_eq!(one.val, 1);
        assert_eq!(one / base * (base * base) - base / one, zero);
    }
}
