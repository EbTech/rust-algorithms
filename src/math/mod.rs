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

// TODO: deduplicate modular arithmetic code with num::Field
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
fn mod_exp(mut base: i64, mut exp: u64, m: i64) -> i64 {
    assert!(m >= 1);
    let mut ans = 1 % m;
    base %= m;
    while exp > 0 {
        if exp % 2 == 1 {
            ans = mod_mul(ans, base, m);
        }
        base = mod_mul(base, base, m);
        exp /= 2;
    }
    pos_mod(ans, m)
}

fn is_strong_probable_prime(n: i64, exp: u64, r: i64, a: i64) -> bool {
    let mut x = mod_exp(a, exp, n);
    if x == 1 || x == n - 1 {
        return true;
    }
    for _ in 1..r {
        x = mod_mul(x, x, n);
        if x == n - 1 {
            return true;
        }
    }
    false
}

/// Assuming x >= 0, returns whether x is prime
pub fn is_prime(n: i64) -> bool {
    const BASES: [i64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    assert!(n >= 0);
    match n {
        0 | 1 => false,
        2 | 3 => true,
        _ if n % 2 == 0 => false,
        _ => {
            let r = (n - 1).trailing_zeros() as i64;
            let exp = (n - 1) as u64 >> r;
            BASES
                .iter()
                .all(|&base| base > n - 2 || is_strong_probable_prime(n, exp, r, base))
        }
    }
}

fn pollard_rho(n: i64) -> i64 {
    for a in 1..n {
        let f = |x| pos_mod(mod_mul(x, x, n) + a, n);
        let mut x = 2;
        let mut y = 2;
        loop {
            x = f(x);
            y = f(f(y));
            let div = num::fast_gcd(x - y, n);
            if div == n {
                break;
            } else if div > 1 {
                return div;
            }
        }
    }
    panic!("No divisor found!");
}

/// Assuming x >= 1, finds the prime factorization of n
/// TODO: pollard_rho needs randomization to ensure correctness in contest settings!
pub fn factorize(n: i64) -> Vec<i64> {
    assert!(n >= 1);
    let r = n.trailing_zeros() as usize;
    let mut factors = vec![2; r];
    let mut stack = match n >> r {
        1 => vec![],
        x => vec![x],
    };
    while let Some(top) = stack.pop() {
        if is_prime(top) {
            factors.push(top);
        } else {
            let div = pollard_rho(top);
            stack.push(div);
            stack.push(top / div);
        }
    }
    factors.sort_unstable();
    factors
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
    fn test_miller() {
        assert_eq!(is_prime(2), true);
        assert_eq!(is_prime(4), false);
        assert_eq!(is_prime(6), false);
        assert_eq!(is_prime(8), false);
        assert_eq!(is_prime(269), true);
        assert_eq!(is_prime(1000), false);
        assert_eq!(is_prime(1_000_000_007), true);
        assert_eq!(is_prime((1 << 61) - 1), true);
        assert_eq!(is_prime(7156857700403137441), false);
    }

    #[test]
    fn test_pollard() {
        assert_eq!(factorize(1), vec![]);
        assert_eq!(factorize(2), vec![2]);
        assert_eq!(factorize(4), vec![2, 2]);
        assert_eq!(factorize(12), vec![2, 2, 3]);
        assert_eq!(
            factorize(7156857700403137441),
            vec![11, 13, 17, 19, 29, 37, 41, 43, 61, 97, 109, 127]
        );
    }
}
