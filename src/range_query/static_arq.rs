//! Associative Range Query Tree based on [Al.Cash's compact representation]
//! (http://codeforces.com/blog/entry/18051).
use super::ArqSpec;

/// Colloquially known as a "segtree" in the sport programming literature, it
/// represents a sequence of elements a_i (0 <= i < size) from a monoid (M, +)
/// on which we want to support fast range operations:
///
/// - modify(l, r, f) replaces a_i (l <= i <= r) by f(a_i) for an endomorphism f
/// - query(l, r) returns the aggregate a_l + a_{l+1} + ... + a_r
///
/// Future work: ArqTree would lend itself naturally to Rust's ownership system.
/// Initially, we should only have access to the root nodes:
///            if size is a power of two, there is a unique root at index 1.
/// arq.push(i) locks i and acquires access to its children.
/// arq.pull(i) is called when the lock on i is released.
pub struct StaticArq<T: ArqSpec> {
    val: Vec<T::M>,
    app: Vec<Option<T::F>>,
}

impl<T: ArqSpec> StaticArq<T> {
    /// Initializes a static balanced tree on top of the given sequence.
    pub fn new(init_val: Vec<T::M>) -> Self {
        let size = init_val.len();
        let mut val = (0..size).map(|_| T::identity()).collect::<Vec<_>>();
        val.append(&mut { init_val });
        let app = vec![None; size];

        let mut arq = Self { val, app };
        for p in (0..size).rev() {
            arq.pull(p);
        }
        arq
    }

    fn apply(&mut self, p: usize, f: &T::F) {
        self.val[p] = T::apply(f, &self.val[p]);
        if let Some(lazy) = self.app.get_mut(p) {
            let h = match *lazy {
                Some(ref g) => T::compose(f, g),
                None => f.clone(),
            };
            *lazy = Some(h);
        }
    }

    fn push(&mut self, p: usize) {
        if let Some(ref f) = self.app[p].take() {
            self.apply(p << 1, f);
            self.apply(p << 1 | 1, f);
        }
    }

    fn pull(&mut self, p: usize) {
        self.val[p] = T::op(&self.val[p << 1], &self.val[p << 1 | 1]);
    }

    fn push_to(&mut self, p: usize) {
        let one_plus_floor_log_p = (p + 1).next_power_of_two().trailing_zeros();
        for s in (1..one_plus_floor_log_p).rev() {
            self.push(p >> s);
        }
    }

    fn pull_from(&mut self, mut p: usize) {
        while p > 1 {
            p >>= 1;
            self.pull(p);
        }
    }

    /// Applies the endomorphism f to all entries from l to r, inclusive.
    /// If l == r, the updates are eager. Otherwise, they are lazy.
    ///
    /// # Panics
    ///
    /// Panics if r >= size. Note that l > r is valid, meaning an empty range.
    pub fn modify(&mut self, mut l: usize, mut r: usize, f: &T::F) {
        l += self.app.len();
        r += self.app.len();
        if l < r {
            self.push_to(l);
        }
        self.push_to(r);
        let (mut l0, mut r0) = (1, 1);
        while l <= r {
            if l & 1 == 1 {
                self.apply(l, f);
                l0 = l0.max(l);
                l += 1;
            }
            if r & 1 == 0 {
                self.apply(r, f);
                r0 = r0.max(r);
                r -= 1;
            }
            l >>= 1;
            r >>= 1;
        }
        self.pull_from(l0);
        self.pull_from(r0);
    }

    /// Returns the aggregate range query on all entries from l to r, inclusive.
    ///
    /// # Panics
    ///
    /// Panics if r >= size. Note that l > r is valid, meaning an empty range.
    pub fn query(&mut self, mut l: usize, mut r: usize) -> T::M {
        l += self.app.len();
        r += self.app.len();
        if l < r {
            self.push_to(l);
        }
        self.push_to(r);
        let (mut l_agg, mut r_agg) = (T::identity(), T::identity());
        while l <= r {
            if l & 1 == 1 {
                l_agg = T::op(&l_agg, &self.val[l]);
                l += 1;
            }
            if r & 1 == 0 {
                r_agg = T::op(&self.val[r], &r_agg);
                r -= 1;
            }
            l >>= 1;
            r >>= 1;
        }
        T::op(&l_agg, &r_agg)
    }
}

/// An example of binary search on the tree of a StaticArq.
/// In this case, we use RMQ to locate the leftmost negative element.
/// To ensure the existence of a valid root note (i == 1) from which to descend,
/// the tree's size must be a power of two.
pub fn first_negative(arq: &mut StaticArq<super::specs::AssignMin>) -> Option<usize> {
    assert!(arq.app.len().is_power_of_two());
    let mut p = 1;
    if arq.val[p] >= 0 {
        None
    } else {
        while p < arq.app.len() {
            arq.push(p);
            p <<= 1;
            if arq.val[p] >= 0 {
                p |= 1;
            }
        }
        Some(p - arq.app.len())
    }
}
