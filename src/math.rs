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

/// Finds (d, x, y) such that d = gcd(a, b) = ax + by.
pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a.abs(), a.signum(), 0)
    } else {
        let (d, x, y) = extended_gcd(b, a % b);
        (d, y, x - y * (a / b))
    }
}

/// Assuming a != 0, finds smallest y >= 0 such that ax + by = c.
///
/// # Panics
///
/// Panics if a == 0.
pub fn canon_egcd(a: i64, b: i64, c: i64) -> Option<(i64, i64, i64)> {
    let (d, _, yy) = extended_gcd(a, b);
    if c % d == 0 {
        let z = (a / d).abs();
        let y = (yy * (c / d) % z + z) % z;
        let x = (c - b * y) / a;
        Some((d, x, y))
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
    }
}
