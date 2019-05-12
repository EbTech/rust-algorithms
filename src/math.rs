//! Number-theoretic utilities for contest problems.

/// Modular exponentiation by repeated squaring: returns base^exp % m.
///
/// # Panics
///
/// Panics if m == 0. May panic on overflow if m * m > 2^63.
pub fn mod_pow(mut base: u64, mut exp: u64, m: u64) -> u64 {
    let mut result = 1 % m;
    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % m;
        }
        base = (base * base) % m;
        exp /= 2;
    }
    result
}

/// Finds (d, coef_a, coef_b) such that d = gcd(a, b) = a * coef_a + b * coef_b.
pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a.abs(), a.signum(), 0)
    } else {
        let (d, coef_a, coef_b) = extended_gcd(b, a % b);
        (d, coef_b, coef_a - coef_b * (a / b))
    }
}

pub fn extended_gcd_by_iter(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a.abs(), a.signum(), 0)
    } else {
        let mut c;
        let mut temp;
        let (mut d, mut q, mut t) = (b, a / b, a % b);
        let (mut x0, mut y0, mut x1, mut y1) = (1_i64, 0_i64, 0_i64, 1_i64);

        while t != 0 {
            temp = x0;
            x0 = x1;
            x1 = temp - q * x1;
            temp = y0;
            y0 = y1;
            y1 = temp - q * y1;
            c = d;
            d = t;
            q = c / d;
            t = c % d;
        }
        (d, x1, y1)
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
    fn test_mod_inverse() {
        let p = 1_000_000_007;
        let base = 31;

        let base_inv = mod_pow(base, p - 2, p);
        let identity = (base * base_inv) % p;

        assert_eq!(identity, 1);
    }

    #[test]
    fn test_egcd() {
        let (a, b) = (14, 35);

        let (d, x, y) = extended_gcd(a, b);
        assert_eq!(d, 7);
        assert_eq!(a * x + b * y, d);

        assert_eq!(canon_egcd(a, b, d), Some((d, -2, 1)));
        assert_eq!(canon_egcd(b, a, d), Some((d, -1, 3)));

        let (d, x, y) = extended_gcd_by_iter(a, b);
        assert_eq!(d, 7);
        assert_eq!(a * x + b * y, d);

        assert_eq!(canon_egcd(a, b, d), Some((d, -2, 1)));
        assert_eq!(canon_egcd(b, a, d), Some((d, -1, 3)));
    }
}
