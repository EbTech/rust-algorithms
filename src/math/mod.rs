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

fn pos_mod(n: i64, m: i64) -> i64 {
    if n < 0 {
        n + m
    } else {
        n
    }
}

fn mod_mul(a: i64, b: i64, m: i64) -> i64 {
    pos_mod((a as i128 * b as i128 % m as i128) as i64, m)
}

/// Assuming m >= 2 and exp >= 0, finds base ^ exp % m in logarithmic time
pub fn mod_exp(mut base: i64, mut exp: i64, m: i64) -> i64 {
    assert!(m >= 2);
    assert!(exp >= 0);
    let mut ans = 1 % m;
    base = base % m;
    while exp > 0 {
        if exp % 2 == 1 {
            ans = mod_mul(ans, base, m);
        }
        base = mod_mul(base, base, m);
        exp /= 2;
    }
    pos_mod(ans, m)
}

/// Assuming m >= 2, finds multiplicative inverse of n under modulus m
pub fn mod_inv(n: i64, m: i64) -> i64 {
    mod_exp(n, m - 2, m)
}

fn miller_test(n: i64, d: i64, r: i64, a: i64) -> bool {
    let mut x = mod_exp(a, d, n);
    if x == 1 || x == n - 1 {
        return true;
    }
    for _ in 0..r {
        x = mod_mul(x, x, n);
        if x == n - 1 {
            return true;
        }
    }
    false
}

const BASES: [i64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
/// Assuming x >= 0, returns true if x is prime
pub fn is_prime(n: i64) -> bool {
    assert!(n >= 0);
    if n <= 1 {
        return false;
    }
    if n <= 3 {
        return true;
    }
    let mut d = n - 1;
    let mut r = 0;
    while d % 2 == 0 {
        d /= 2;
        r += 1;
    }
    for base in BASES.iter() {
        if *base <= n - 2 && !miller_test(n, d, r, *base) {
            return false;
        }
    }
    true
}

fn pollard_rho(n: i64) -> i64 {
    for a in 1..n {
        let f = |x| pos_mod(mod_mul(x, x, n) + a, n);
        for b in 0..n {
            let mut x = b;
            let mut y = b;
            loop {
                x = f(x);
                y = f(f(y));
                let p = num::fast_gcd((x - y).abs(), n);
                if p > 1 && p < n {
                    return p;
                }
                if x == y {
                    break;
                }
            }
        }
    }
    panic!("No divisor found!");
}

pub fn factorize(n: i64) -> Vec<i64> {
    Vec::new()
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

    #[test]
    fn test_modexp() {
        let m = 1_000_000_007;
        assert_eq!(mod_exp(0, 0, m), 1);
        assert_eq!(mod_exp(0, 1, m), 0);
        assert_eq!(mod_exp(0, 10, m), 0);
        assert_eq!(mod_exp(123, 456, m), 565291922);
    }

    #[test]
    fn test_modinv() {
        let m = 1_000_000_007;
        assert_eq!(mod_inv(1, m), 1);
        assert_eq!(mod_inv(-1, m), m - 1);
        assert_eq!(mod_inv(3301, m), 756740387);
        assert_eq!(mod_inv(756740387, m), 3301);
    }

    #[test]
    fn test_miller() {
        assert_eq!(is_prime(2), true);
        assert_eq!(is_prime(4), false);
        assert_eq!(is_prime(269), true);
        assert_eq!(is_prime(1_000_000_007), true);
        assert_eq!(is_prime((1 << 61) - 1), true);
        assert_eq!(is_prime(7156857700403137441), false);
    }

    #[test]
    fn test_pollard() {
        println!("{}", pollard_rho(4));
    }
}
