

// Find x,y such that d = gcd(a,b) = ax + by
// * a^-1 (mod m): let (d, x, y) = egcd(a,m), 1); assert!(d, 1); return (x+m)%m;
pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a.abs(), a.signum(), 0)
    }
    else {
        let (d, x, y) = extended_gcd(b, a % b);
        (d, y, x - y * (a / b))
    }
}

// Assuming a != 0, find smallest y >= 0 such that ax + by = c (if possible)
pub fn canon_egcd(a: i64, b: i64, c: i64) -> Option<(i64, i64, i64)> {
    let (d, _, mut y) = extended_gcd(a, b);
    let z = (a / d).abs();
    if c % d != 0 {
        None
    }
    else {
        y = (y*(c/d)%z + z)%z;
        let x = (c - b*y)/a;
        Some((d, x, y))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test_gcd() {
        let (a, b) = (14, 35);
        let (d, x, y) = extended_gcd(a, b);
        assert_eq!(d, 7);
        assert_eq!(a*x + b*y, d);
        assert_eq!(canon_egcd(a, b, d), Some((d, -2, 1)));
        assert_eq!(canon_egcd(b, a, d), Some((d, -1, 3)));
    }
}
