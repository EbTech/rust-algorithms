//! Associative Range Query Tree based on [Al.Cash's compact representation]
//! (http://codeforces.com/blog/entry/18051).

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
pub struct ArqTree<T: ArqSpec> {
    app: Vec<Option<T::F>>,
    val: Vec<T::M>,
}

impl<T: ArqSpec> ArqTree<T> {
    /// Initializes a static balanced tree on top of the given sequence.
    pub fn new(init_val: Vec<T::M>) -> Self {
        let size = init_val.len();
        let mut val = (0..size).map(|_| T::identity()).collect::<Vec<_>>();
        val.append(&mut { init_val });
        let app = vec![None; size];

        let mut arq = Self { app, val };
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

pub trait ArqSpec {
    /// Type of data representing an endomorphism.
    // Note that while a Fn(M) -> M may seem like a more natural representation
    // for an endomorphism, compositions would then have to delegate to each of
    // their parts. This representation is more efficient.
    type F: Clone;
    /// Type of monoid elements.
    type M;

    /// For eager updates, compose() ho be unimplemented!(). For lazy updates:
    /// Require for all f,g,a: apply(compose(f, g), a) = apply(f, apply(g, a))
    fn compose(f: &Self::F, g: &Self::F) -> Self::F;
    /// For eager updates, apply() can assume to act on a leaf. For lazy updates:
    /// Require for all f,a,b: apply(f, op(a, b)) = op(apply(f, a), apply(f, b))
    fn apply(f: &Self::F, a: &Self::M) -> Self::M;
    /// Require for all a,b,c: op(a, op(b, c)) = op(op(a, b), c)
    fn op(a: &Self::M, b: &Self::M) -> Self::M;
    /// Require for all a: op(a, identity()) = op(identity(), a) = a
    fn identity() -> Self::M;
}

/// Range Minimum Query (RMQ), a classic application of ARQ.
/// modify(l, r, &f) sets all entries a[l..=r] to f.
/// query(l, r) finds the minimum value in a[l..=r].
//
// Exercises: try augmenting this struct to find the index of a minimum element
// in a range query, as well as the number of elements equal to the minimum.
// Then instead of overwriting values with a constant assignment a[i] = f,
// try supporting addition: a[i] += f.
pub enum AssignMin {}
impl ArqSpec for AssignMin {
    type F = i64;
    type M = i64;
    fn compose(&f: &Self::F, _: &Self::F) -> Self::F {
        f
    }
    fn apply(&f: &Self::F, _: &Self::M) -> Self::M {
        f
    }
    fn op(&a: &Self::M, &b: &Self::M) -> Self::M {
        a.min(b)
    }
    fn identity() -> Self::M {
        Self::M::max_value()
    }
}

/// An example of binary search on an ArqTree.
/// In this case, we use RMQ to locate the leftmost negative element.
/// To ensure the existence of a valid root note (i == 1) from which to descend,
/// the tree size must be a power of two.
pub fn first_negative(arq: &mut ArqTree<AssignMin>) -> i32 {
    assert!(arq.app.len().is_power_of_two());
    let mut i = 1;
    if arq.val[i] >= 0 {
        return -1;
    }
    while i < arq.app.len() {
        arq.push(i);
        i <<= 1;
        if arq.val[i] >= 0 {
            i |= 1;
        }
    }
    let pos = i - arq.app.len();
    pos as i32
}

/// Range Sum Query, a slightly trickier classic application of ARQ.
/// modify(l, r, &f) sets all entries a[l..=r] to f.
/// query(l, r) sums all the entries a[l..=r].
///
/// # Panics
///
/// Associated functions will panic on overflow.
//
// Note that the apply() operation on raw entries is undefined: while leaf nodes
// should simply be set to f, internal nodes must be set to f * size_of_subtree.
// Thus, our monoid type M should store the pair (entry, size_of_subtree).
//
// In mathematical jargon, we say that constant assignment f(a) = f is not an
// endomorphism on (i64, +) because f(a+b) = f != 2*f = f(a) + f(b).
// On the other hand, f((a, s)) = (f*s, s) is indeed an endomorphism on pairs
// with vector addition: f((a, s) + (b, t)) = f((a+b, s+t)) = (f*(s+t), s+t)
//                       = (f*s, s) + (f*t, t) = f((a,s)) + f((b,t)).
pub enum AssignSum {}
impl ArqSpec for AssignSum {
    type F = i64;
    type M = (i64, i64);
    fn compose(&f: &Self::F, _: &Self::F) -> Self::F {
        f
    }
    fn apply(&f: &Self::F, &(_, s): &Self::M) -> Self::M {
        (f * s, s)
    }
    fn op(&(a, s): &Self::M, &(b, t): &Self::M) -> Self::M {
        (a + b, s + t)
    }
    fn identity() -> Self::M {
        (0, 0)
    }
}

/// Supply & Demand, based on https://codeforces.com/gym/102218/problem/F
/// modify(i, i, &(p, o)) increases supply by p and demand by o at time i.
/// query(l, r) computes total supply and demand at times l to r, as well as
//              how much of the supply is subsequently met by the demand.
//
// Note that the apply() operation is only correct when applied to leaf nodes.
// Therefore, modify() must only be used in "eager" mode, i.e., with l == r.
// compose() should be unimplemented!() to prevent accidental "lazy" updates.
pub enum SupplyDemand {}
impl ArqSpec for SupplyDemand {
    type F = (i64, i64);
    type M = (i64, i64, i64); // production, orders, sales
    fn compose(_: &Self::F, _: &Self::F) -> Self::F {
        unimplemented!()
    }
    fn apply(&(p_add, o_add): &Self::F, &(p, o, _): &Self::M) -> Self::M {
        let p = p + p_add;
        let o = o + o_add;
        (p, o, p.min(o))
    }
    fn op((p1, o1, s1): &Self::M, (p2, o2, s2): &Self::M) -> Self::M {
        let extra = (p1 - s1).min(o2 - s2);
        (p1 + p2, o1 + o2, s1 + s2 + extra)
    }
    fn identity() -> Self::M {
        (0, 0, 0)
    }
}

/// A generic implementation of Mo's algorithm, aka Query Sqrt Decomposition.
/// It answers q offline queries over intervals in 0..n by shifting the query
/// interval's endpoints by one position at a time.
/// Each endpoint incurs a total cost of at most n * sqrt(q * L_OP * R_OP).
pub trait MoState {
    type Q;
    type A;

    /// cost ratio L_OP / R_OP between a left endpoint and a right endpoint move
    const L_R_RATIO: f64 = 1.0;

    fn query(&self, q: &Self::Q) -> Self::A;
    fn insert_left(&mut self, pos: usize);
    fn remove_left(&mut self, pos: usize);

    fn insert_right(&mut self, pos: usize) {
        self.insert_left(pos);
    }
    fn remove_right(&mut self, pos: usize) {
        self.remove_left(pos);
    }
    /// After initializing self to a state corresponding to an empty interval,
    /// call this function to answer all your queries.
    fn process(&mut self, queries: &[(usize, usize, Self::Q)]) -> Vec<Self::A> {
        let q = queries.len();
        let mut q_positions: Vec<usize> = (0..q).collect();
        if let Some(max_r) = queries.iter().map(|&(_, r, _)| r).max() {
            let q_adjusted = q as f64 * Self::L_R_RATIO;
            let bucket_width = 1 + max_r / q_adjusted.sqrt() as usize;
            q_positions.sort_unstable_by_key(|&i| {
                let (l, mut r) = (queries[i].0, queries[i].1);
                let bucket = l / bucket_width;
                if bucket % 2 != 0 {
                    r = max_r - r;
                }
                (bucket, r)
            });
        }

        let (mut cur_l, mut cur_r) = (1, 0);
        let mut answers = Vec::with_capacity(queries.len());
        for i in q_positions {
            let (l, r, ref q) = queries[i];
            while cur_l > l {
                cur_l -= 1;
                self.insert_left(cur_l);
            }
            while cur_r < r {
                cur_r += 1;
                self.insert_right(cur_r);
            }
            while cur_l < l {
                self.remove_left(cur_l);
                cur_l += 1;
            }
            while cur_r > r {
                self.remove_right(cur_r);
                cur_r -= 1;
            }
            answers.push((i, self.query(q)));
        }
        answers.sort_unstable_by_key(|&(i, _)| i);
        answers.into_iter().map(|(_, ans)| ans).collect()
    }
}

pub struct DistinctVals {
    vals: Vec<usize>,
    counts: Vec<usize>,
    distinct: usize,
}
impl DistinctVals {
    pub fn new(vals: Vec<usize>) -> Self {
        let max_val = vals.iter().cloned().max().unwrap_or(0);
        Self {
            vals,
            counts: vec![0; max_val + 1],
            distinct: 0,
        }
    }
}
impl MoState for DistinctVals {
    type Q = ();
    type A = usize;
    fn query(&self, _: &Self::Q) -> Self::A {
        self.distinct
    }
    fn insert_left(&mut self, pos: usize) {
        let v = self.vals[pos];
        if self.counts[v] == 0 {
            self.distinct += 1;
        }
        self.counts[v] += 1;
    }
    fn remove_left(&mut self, pos: usize) {
        let v = self.vals[pos];
        self.counts[v] -= 1;
        if self.counts[v] == 0 {
            self.distinct -= 1;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_rmq() {
        let mut arq = ArqTree::<AssignMin>::new(vec![0; 10]);

        assert_eq!(arq.query(0, 9), 0);

        arq.modify(2, 4, &-5);
        arq.modify(5, 7, &-3);
        arq.modify(1, 6, &1);

        assert_eq!(arq.query(0, 9), -3);
    }

    #[test]
    fn test_rmq_binary_search() {
        let vec = vec![0, 1, -2, 3, -4, -5, 6, -7];
        let mut arq = ArqTree::<AssignMin>::new(vec);
        let pos = first_negative(&mut arq);

        arq.modify(2, 7, &0);
        let pos_zeros = first_negative(&mut arq);

        assert_eq!(pos, 2);
        assert_eq!(pos_zeros, -1);
    }

    #[test]
    fn test_range_sum() {
        let mut arq = ArqTree::<AssignSum>::new(vec![(0, 1); 10]);

        assert_eq!(arq.query(0, 9), (0, 10));

        arq.modify(1, 3, &10);
        arq.modify(3, 5, &1);

        assert_eq!(arq.query(0, 9), (23, 10));
        assert_eq!(arq.query(10, 4), (0, 0));
    }

    #[test]
    fn test_supply_demand() {
        let mut arq = ArqTree::<SupplyDemand>::new(vec![(0, 0, 0); 10]);

        arq.modify(1, 1, &(25, 100));
        arq.modify(3, 3, &(100, 30));
        arq.modify(9, 9, &(0, 20));

        assert_eq!(arq.query(0, 9), (125, 150, 75));
    }

    #[test]
    fn test_mos_algorithm() {
        let queries = vec![(0, 2, ()), (5, 5, ()), (2, 6, ()), (0, 6, ())];
        let arr = vec![4, 8, 4, 7, 1, 9, 8];

        let answers = DistinctVals::new(arr).process(&queries);

        assert_eq!(answers, vec![2, 1, 5, 5]);
    }
}
