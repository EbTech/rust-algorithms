//! Associative Range Query Tree based on [Al.Cash's compact representation]
//! (http://codeforces.com/blog/entry/18051).

/// Colloquially known as a "segtree" in the sport programming literature, it
/// represents a sequence of elements a_i (0 <= i < size) from a monoid (M, +) on
/// which we want to support fast range operations:
///
/// - modify(l, r, f) replaces a_i (l <= i <= r) by f(a_i) for a homomorphism f
/// - query(l, r) returns the aggregate a_l + a_{l+1} + ... + a_r
///
/// To customize, simply change the commented lines.
/// In this example, we chose to support range sum queries and range constant
/// assignments. Since constant assignment f_c(a) = c is not a homomorphism over
/// integers, we have to augment the monoid type, using the 2D vector (a_i, 1)
/// instead of a_i. You may check that f_c((a, s)) = (c*s, s) is a homomorphism.
pub struct ArqTree {
    d: Vec<Option<i64>>,
    t: Vec<i64>,
    s: Vec<i64>,
}

impl ArqTree {
    /// Initializes a sequence of identity elements.
    pub fn new(size: usize) -> Self {
        let mut s = vec![1; 2 * size];
        for i in (0..size).rev() {
            s[i] = s[i << 1] + s[i << 1 | 1];
        }
        Self {
            d: vec![None; size],
            t: vec![0; 2 * size], // monoid identity
            s: s,
        }
    }

    fn apply(&mut self, p: usize, f: i64) {
        self.t[p] = f * self.s[p]; // hom application
        if p < self.d.len() {
            self.d[p] = Some(f); // hom composition
        }
    }

    fn push(&mut self, p: usize) {
        for s in (1..32).rev() {
            let i = p >> s;
            if let Some(f) = self.d[i] {
                self.apply(i << 1, f);
                self.apply(i << 1 | 1, f);
                self.d[i] = None;
            }
        }
    }

    fn pull(&mut self, mut p: usize) {
        while p > 1 {
            p >>= 1;
            if self.d[p] == None {
                self.t[p] = self.t[p << 1] + self.t[p << 1 | 1]; // monoid op
            }
        }
    }

    /// Performs the homomorphism f on all entries from l to r, inclusive.
    pub fn modify(&mut self, mut l: usize, mut r: usize, f: i64) {
        l += self.d.len();
        r += self.d.len();
        let (l0, r0) = (l, r);
        self.push(l0);
        self.push(r0);
        while l <= r {
            if l & 1 == 1 {
                self.apply(l, f);
                l += 1;
            }
            if r & 1 == 0 {
                self.apply(r, f);
                r -= 1;
            }
            l >>= 1;
            r >>= 1;
        }
        self.pull(l0);
        self.pull(r0);
    }

    /// Returns the aggregate range query on all entries from l to r, inclusive.
    pub fn query(&mut self, mut l: usize, mut r: usize) -> i64 {
        l += self.d.len();
        r += self.d.len();
        self.push(l);
        self.push(r);
        let mut res = 0; // monoid identity
        while l <= r {
            if l & 1 == 1 {
                res = res + self.t[l]; // monoid op
                l += 1;
            }
            if r & 1 == 0 {
                res = self.t[r] + res; // monoid op
                r -= 1;
            }
            l >>= 1;
            r >>= 1;
        }
        res
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_arq_tree() {
        let mut arq = ArqTree::new(10);

        arq.modify(1, 3, 10);
        arq.modify(3, 5, 1);

        assert_eq!(arq.query(0, 9), 23);
    }
}
