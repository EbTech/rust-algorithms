//! Number-theoretic utilities for contest problems.
pub mod fft;
pub mod num;

/// Finds (d, coef_a, coef_b) such that d = gcd(a, b) = a * coef_a + b * coef_b.
pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a.abs(), a.signum(), 0)
    } else {
        let (d, coef_b, coef_a) = extended_gcd(b, a % b);
        (d, coef_a, coef_b - coef_a * (a / b))
    }
}

/// Assuming a != 0, finds smallest coef_b >= 0 such that a * coef_a + b * coef_b = c.
///
/// # Panics
///
/// Panics if a == 0.
pub fn canon_egcd(a: i64, b: i64, c: i64) -> Option<(i64, i64, i64)> {
    let (d, _, coef_b_init) = extended_gcd(a, b);
    if c % d == 0 {
        let a_d = (a / d).abs();
        let coef_b = (coef_b_init * (c / d) % a_d + a_d) % a_d;
        let coef_a = (c - b * coef_b) / a;
        Some((d, coef_a, coef_b))
    } else {
        None
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_egcd() {
        let (a, b) = (14, 35);

        let (d, x, y) = extended_gcd(a, b);
        assert_eq!(d, 7);
        assert_eq!(a * x + b * y, d);

        assert_eq!(canon_egcd(a, b, d), Some((d, -2, 1)));
        assert_eq!(canon_egcd(b, a, d), Some((d, -1, 3)));
    }
}
