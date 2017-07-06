//! Associative Range Query Tree based on [Al.Cash's compact representation]
//! (http://codeforces.com/blog/entry/18051).

/// Colloquially known as a "segtree" in the sport programming literature, it
/// represents a sequence of elements a_i (0 <= i < size) from a monoid (M, +)
/// on which we want to support fast range operations:
///
/// - modify(l, r, f) replaces a_i (l <= i <= r) by f(a_i) for an endomorphism f
/// - query(l, r) returns the aggregate a_l + a_{l+1} + ... + a_r
pub struct ArqTree<T: ArqSpec> {
    d: Vec<Option<T::F>>,
    t: Vec<T::M>,
}

impl<T: ArqSpec> ArqTree<T>
where
    T::F: Clone,
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
            d: vec![None; size],
            t: t,
        }
    }

    fn apply(&mut self, p: usize, f: &T::F) {
        self.t[p] = T::apply(f, &self.t[p]);
        if p < self.d.len() {
            let h = match self.d[p] {
                Some(ref g) => T::compose(f, g),
                None => f.clone(),
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

    /// Applies the endomorphism f to all entries from l to r, inclusive.
    ///
    /// # Panics
    ///
    /// Panics if l or r is out of range.
    pub fn modify(&mut self, mut l: usize, mut r: usize, f: &T::F) {
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
    ///
    /// # Panics
    ///
    /// Panics if l or r is out of range.
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
    /// Type of data representing an endomorphism.
    // Note that while a Fn(M) -> M may seem like a more natural representation
    // for an endomorphism, compositions would then have to delegate to each of
    // their parts. This representation is more efficient.
    type F;
    /// Type of monoid elements.
    type M;
    /// Require for all f,g,a: apply(compose(f, g), a) = apply(f, apply(g, a))
    fn compose(f: &Self::F, g: &Self::F) -> Self::F;
    /// Require for all f,a,b: apply(f, op(a, b)) = op(apply(f, a), apply(f, b))
    fn apply(f: &Self::F, a: &Self::M) -> Self::M;
    /// Require for all a,b,c: op(a, op(b, c)) = op(op(a, b), c)
    fn op(a: &Self::M, b: &Self::M) -> Self::M;
    /// Require for all a: op(a, identity()) = op(identity(), a) = a
    fn identity() -> Self::M;
}

/// In this example, we want to support range sum queries and range constant
/// assignments. Note that constant assignment f_c(a) = c is not a endomorphism
/// on (i64, +) because f_c(a+b) = c != 2*c = f_c(a) + f_c(b). In intuitive
/// terms, the problem is that the internal nodes of the tree should really be
/// set to a multiple of c, corresponding to the subtree size. So let's augment
/// the monoid type with size information, using the 2D vector (a_i,1) instead
/// of a_i. Now check that f_c((a, s)) = (c*s, s) is indeed an endomorphism on
/// vector addition: f_c((a,s)+(b,t)) = f_c((a+b,s+t)) = (c*(s+t),s+t)
///                 = (c*s,s)+(c*t,t) = f_c((a,s)) + f_c((b,t)).
///
/// # Panics
///
/// Associated functions will panic on overflow.
pub struct AssignSum;

impl ArqSpec for AssignSum {
    type F = i64;
    type M = (i64, i64);
    fn compose(f: &Self::F, _: &Self::F) -> Self::F {
        *f
    }
    fn apply(f: &Self::F, a: &Self::M) -> Self::M {
        (f * a.1, a.1)
    }
    fn op(a: &Self::M, b: &Self::M) -> Self::M {
        (a.0 + b.0, a.1 + b.1)
    }
    fn identity() -> Self::M {
        (0, 0)
    }
}

/// Range Minimum Query, a classic form of associative range query.
// Exercises: try augmenting this struct to find the index of a minimum element
// in a range query, as well as the number of elements that match the minimum.
// Then instead of assigning to a range, try to support the operation of
// incrementing each element in a range by a given offset!
pub struct AssignMin;

impl ArqSpec for AssignMin {
    type F = i64;
    type M = i64;
    fn compose(f: &Self::F, _: &Self::F) -> Self::F {
        *f
    }
    fn apply(f: &Self::F, _: &Self::M) -> Self::M {
        *f
    }
    fn op(a: &Self::M, b: &Self::M) -> Self::M {
        ::std::cmp::min(*a, *b)
    }
    fn identity() -> Self::M {
        Self::M::max_value()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_range_sum() {
        let mut arq = ArqTree::<AssignSum>::new(vec![(0, 1); 10]);

        assert_eq!(arq.query(0, 9), (0, 10));

        arq.modify(1, 3, &10);
        arq.modify(3, 5, &1);

        assert_eq!(arq.query(0, 9), (23, 10));
    }

    #[test]
    fn test_rmq() {
        let mut arq = ArqTree::<AssignMin>::new(vec![0; 10]);

        assert_eq!(arq.query(0, 9), 0);

        arq.modify(2, 4, &-5);
        arq.modify(5, 7, &-3);
        arq.modify(1, 6, &1);

        assert_eq!(arq.query(0, 9), -3);
    }
}
