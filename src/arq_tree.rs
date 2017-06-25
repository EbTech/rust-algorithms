//! Associative Range Query Tree based on [Al.Cash's compact representation]
//! (http://codeforces.com/blog/entry/18051).

/// Colloquially known as a "segtree" in the sport programming literature, it
/// represents a sequence of elements a_i (0 <= i < size) from a monoid (M, +) on
/// which we want to support fast range operations:
///
/// - modify(l, r, f) replaces a_i (l <= i <= r) by f(a_i) for a homomorphism f
/// - query(l, r) returns the aggregate a_l + a_{l+1} + ... + a_r
pub struct ArqTree<T: ArqSpec> {
    d: Vec<Option<T::H>>,
    t: Vec<T::M>,
}

impl<T: ArqSpec> ArqTree<T>
where
    T::H: Clone,
{
    /// Initializes a static balanced tree on top of the given sequence.
    pub fn new(mut init: Vec<T::M>) -> Self {
        let size = init.len();
        let mut t = (0..size).map(|_| T::identity()).collect::<Vec<_>>();
        t.append(&mut init);
        for i in (0..size).rev() {
            t[i] = T::op(&t[i << 1], &t[i << 1 | 1]);
        }
        Self {
            d: (0..size).map(|_| None).collect(),
            t: t,
        }
    }

    fn apply(&mut self, p: usize, f: &T::H) {
        self.t[p] = T::apply(f, &self.t[p]);
        if p < self.d.len() {
            let h = if let Some(ref g) = self.d[p] {
                T::compose(f, g)
            } else {
                f.clone()
            };
            self.d[p] = Some(h);
        }
    }

    fn push(&mut self, p: usize) {
        for s in (1..32).rev() {
            let i = p >> s;
            if let Some(ref f) = self.d[i].take() {
                self.apply(i << 1, f);
                self.apply(i << 1 | 1, f);
            }
        }
    }

    fn pull(&mut self, mut p: usize) {
        while p > 1 {
            p >>= 1;
            if self.d[p].is_none() {
                self.t[p] = T::op(&self.t[p << 1], &self.t[p << 1 | 1]);
            }
        }
    }

    /// Performs the homomorphism f on all entries from l to r, inclusive.
    pub fn modify(&mut self, mut l: usize, mut r: usize, f: &T::H) {
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
    pub fn query(&mut self, mut l: usize, mut r: usize) -> T::M {
        l += self.d.len();
        r += self.d.len();
        self.push(l);
        self.push(r);
        let (mut l_agg, mut r_agg) = (T::identity(), T::identity());
        while l <= r {
            if l & 1 == 1 {
                l_agg = T::op(&l_agg, &self.t[l]);
                l += 1;
            }
            if r & 1 == 0 {
                r_agg = T::op(&self.t[r], &r_agg);
                r -= 1;
            }
            l >>= 1;
            r >>= 1;
        }
        T::op(&l_agg, &r_agg)
    }
}

pub trait ArqSpec {
    type H;
    type M;
    fn compose(f: &Self::H, g: &Self::H) -> Self::H;
    fn apply(f: &Self::H, a: &Self::M) -> Self::M;
    fn op(a: &Self::M, b: &Self::M) -> Self::M;
    fn identity() -> Self::M;
}

/// In this example, we chose to support range sum queries and range constant
/// assignments. Since constant assignment f_c(a) = c is not a homomorphism over
/// integers, we have to augment the monoid type, using the 2D vector (a_i, 1)
/// instead of a_i. You may check that f_c((a, s)) = (c*s, s) is a homomorphism.
pub struct AssignAdd;

impl ArqSpec for AssignAdd {
    type H = i64;
    type M = (i64, i64);
    fn compose(f: &Self::H, _: &Self::H) -> Self::H {
        *f
    }
    fn apply(f: &Self::H, a: &Self::M) -> Self::M {
        (f * a.1, a.1)
    }
    fn op(a: &Self::M, b: &Self::M) -> Self::M {
        (a.0 + b.0, a.1 + b.1)
    }
    fn identity() -> Self::M {
        (0, 0)
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_arq_tree() {
        let mut arq = ArqTree::<AssignAdd>::new(vec![(0, 1); 10]);

        assert_eq!(arq.query(0, 9), (0, 10));

        arq.modify(1, 3, &10);
        arq.modify(3, 5, &1);

        assert_eq!(arq.query(0, 9), (23, 10));
    }
}
