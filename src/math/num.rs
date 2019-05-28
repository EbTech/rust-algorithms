//! Rational and Complex numbers, minimally viable for contests.
//! For more optimized and full-featured versions, check out the num crate.
// TODO: Once Rust gets const generics, consider implementing a Mod<p> struct
use super::extended_gcd;
pub use std::f64::consts::PI;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Represents a fraction reduced to lowest terms
#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub struct Rational {
    pub num: i64,
    pub den: i64,
}
impl Rational {
    pub fn new(a: i64, b: i64) -> Self {
        let g = extended_gcd(a, b).0 * b.signum();
        Self {
            num: a / g,
            den: b / g,
        }
    }
    pub fn abs(self) -> Self {
        Self {
            num: self.num.abs(),
            den: self.den,
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
        Self::new(self.num * other.den, self.den * other.num)
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
    pub fn new(a: f64, b: f64) -> Self {
        Self { real: a, imag: b }
    }
    pub fn from_polar(r: f64, th: f64) -> Self {
        Self::new(r * th.cos(), r * th.sin())
    }
    pub fn conjugate(self) -> Self {
        Self::new(self.real, -self.imag)
    }
    pub fn abs_square(self) -> f64 {
        self.real * self.real + self.imag * self.imag
    }
    pub fn argument(self) -> f64 {
        self.imag.atan2(self.real)
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
impl Div for Complex {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let real = self.real * other.real + self.imag * other.imag;
        let imag = self.imag * other.real - self.real * other.imag;
        let denom = other.abs_square();
        Self::new(real / denom, imag / denom)
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
}
