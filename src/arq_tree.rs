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

pub struct DynamicArqNode<T: ArqSpec> {
    val: T::M,
    app: Option<T::F>,
    down: Option<(usize, usize)>,
}

// TODO: can this be replaced by a #[derive(Clone)]?
impl<T: ArqSpec> Clone for DynamicArqNode<T> {
    fn clone(&self) -> Self {
        Self {
            val: T::op(&T::identity(), &self.val),
            app: self.app.clone(),
            down: self.down,
        }
    }
}

impl<T: ArqSpec> DynamicArqNode<T> {
    pub fn new(val: T::M) -> Self {
        Self {
            val,
            app: None,
            down: None,
        }
    }

    fn apply(&mut self, f: &T::F, is_leaf: bool) {
        self.val = T::apply(f, &self.val);
        if !is_leaf {
            let h = match self.app {
                Some(ref g) => T::compose(f, g),
                None => f.clone(),
            };
            self.app = Some(h);
        }
    }
}

/// A dynamic, and optionally persistent, associate range query data structure.
pub struct DynamicArq<T: ArqSpec> {
    l_bound: i64,
    r_bound: i64,
    nodes: Vec<DynamicArqNode<T>>,
    is_persistent: bool,
    initializer: Box<dyn Fn(i64, i64) -> T::M>,
}

impl<T: ArqSpec> DynamicArq<T> {
    pub fn new(
        l_bound: i64,
        r_bound: i64,
        is_persistent: bool,
        initializer: Box<dyn Fn(i64, i64) -> T::M>,
    ) -> Self {
        let val = initializer(l_bound, r_bound);
        let nodes = vec![DynamicArqNode::new(val)];
        Self {
            l_bound,
            r_bound,
            nodes,
            is_persistent,
            initializer,
        }
    }

    pub fn new_with_identity(l_bound: i64, r_bound: i64, is_persistent: bool) -> Self {
        let initializer = Box::new(|_, _| T::identity());
        Self::new(l_bound, r_bound, is_persistent, initializer)
    }

    fn push(&mut self, p: usize, l: i64, r: i64) -> (usize, usize) {
        let m = (l + r) / 2;
        let (lp, rp) = match self.nodes[p].down {
            Some(children) => children,
            None => {
                let l_val = (self.initializer)(l, m);
                let r_val = (self.initializer)(m + 1, r);
                self.nodes.push(DynamicArqNode::new(l_val));
                self.nodes.push(DynamicArqNode::new(r_val));
                (self.nodes.len() - 2, self.nodes.len() - 1)
            }
        };
        if let Some(ref f) = self.nodes[p].app.take() {
            self.nodes[lp].apply(f, l == m);
            self.nodes[rp].apply(f, m + 1 == r);
        }
        (lp, rp)
    }

    fn pull(&mut self, p: usize) {
        let (lp, rp) = self.nodes[p].down.unwrap();
        let left_val = &self.nodes[lp].val;
        let right_val = &self.nodes[rp].val;
        self.nodes[p].val = T::op(left_val, right_val);
    }

    fn clone_node(&mut self, p: usize) -> usize {
        if self.is_persistent {
            let node = self.nodes[p].clone();
            self.nodes.push(node);
            self.nodes.len() - 1
        } else {
            p
        }
    }

    fn m_helper(&mut self, p: usize, l: i64, r: i64, f: &T::F, cl: i64, cr: i64) -> usize {
        if r < cl || cr < l {
            p
        } else if l <= cl && cr <= r /* && self.l == self.r forces eager */ {
            let p_clone = self.clone_node(p);
            self.nodes[p_clone].apply(f, l == r);
            p_clone
        } else {
            let (lp, rp) = self.push(p, cl, cr);
            let cm = (cl + cr) / 2;
            let p_clone = self.clone_node(p);
            let lp_clone = self.m_helper(lp, l, r, f, cl, cm);
            let rp_clone = self.m_helper(rp, l, r, f, cm + 1, cr);
            self.nodes[p_clone].down = Some((lp_clone, rp_clone));
            self.pull(p_clone);
            p_clone
        }
    }

    fn q_helper(&mut self, p: usize, l: i64, r: i64, cl: i64, cr: i64) -> T::M {
        if r < cl || cr < l {
            T::identity()
        } else if l <= cl && cr <= r {
            T::op(&T::identity(), &self.nodes[p].val)
        } else {
            let (lp, rp) = self.push(p, cl, cr);
            let cm = (cl + cr) / 2;
            let l_agg = self.q_helper(lp, l, r, cl, cm);
            let r_agg = self.q_helper(rp, l, r, cm + 1, cr);
            T::op(&l_agg, &r_agg)
        }
    }

    /// Applies the endomorphism f to all entries from l to r, inclusive.
    /// If l == r, the updates are eager. Otherwise, they are lazy.
    pub fn modify(&mut self, p: usize, l: i64, r: i64, f: &T::F) -> usize {
        self.m_helper(p, l, r, f, self.l_bound, self.r_bound)
    }

    /// Returns the aggregate range query on all entries from l to r, inclusive.
    pub fn query(&mut self, p: usize, l: i64, r: i64) -> T::M {
        self.q_helper(p, l, r, self.l_bound, self.r_bound)
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

/// An example of binary search on the tree of a StaticArq.
/// In this case, we use RMQ to locate the leftmost negative element.
/// To ensure the existence of a valid root note (i == 1) from which to descend,
/// the tree's size must be a power of two.
pub fn first_negative_static(arq: &mut StaticArq<AssignMin>) -> i32 {
    assert!(arq.app.len().is_power_of_two());
    let mut p = 1;
    if arq.val[p] >= 0 {
        return -1;
    }
    while p < arq.app.len() {
        arq.push(p);
        p <<= 1;
        if arq.val[p] >= 0 {
            p |= 1;
        }
    }
    let pos = p - arq.app.len();
    pos as i32
}

/// An example of binary search on the tree of a DynamicArq.
/// The tree may have any size, not necessarily a power of two.
pub fn first_negative_dynamic(arq: &mut DynamicArq<AssignMin>, p: usize, cl: i64, cr: i64) -> i64 {
    if arq.nodes[p].val >= 0 {
        -1
    } else if cl == cr {
        cl
    } else {
        let (lp, rp) = arq.push(p, cl, cr);
        let cm = (cl + cr) / 2;
        if arq.nodes[lp].val < 0 {
            first_negative_dynamic(arq, lp, cl, cm)
        } else {
            first_negative_dynamic(arq, rp, cm + 1, cr)
        }
    }
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
        let &max_val = vals.iter().max().unwrap_or(&0);
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
        let mut arq = StaticArq::<AssignMin>::new(vec![0; 10]);

        assert_eq!(arq.query(0, 9), 0);

        arq.modify(2, 4, &-5);
        arq.modify(5, 7, &-3);
        arq.modify(1, 6, &1);

        assert_eq!(arq.query(0, 9), -3);
    }

    #[test]
    fn test_dynamic_rmq() {
        let initializer = Box::new(|_, _| 0);
        let mut arq = DynamicArq::<AssignMin>::new(0, 9, false, initializer);

        assert_eq!(arq.query(0, 0, 9), 0);

        arq.modify(0, 2, 4, &-5);
        arq.modify(0, 5, 7, &-3);
        arq.modify(0, 1, 6, &1);

        assert_eq!(arq.query(0, 0, 9), -3);
    }

    #[test]
    fn test_persistent_rmq() {
        let initializer = Box::new(|_, _| 0);
        let mut arq = DynamicArq::<AssignMin>::new(0, 9, true, initializer);

        let mut p = 0;
        p = arq.modify(p, 2, 4, &-5);
        let snapshot = p;
        p = arq.modify(p, 5, 7, &-3);
        p = arq.modify(p, 1, 6, &1);

        assert_eq!(arq.query(0, 0, 9), 0);
        assert_eq!(arq.query(snapshot, 0, 9), -5);
        assert_eq!(arq.query(p, 0, 9), -3);
    }

    #[test]
    fn test_range_sum() {
        let mut arq = StaticArq::<AssignSum>::new(vec![(0, 1); 10]);

        assert_eq!(arq.query(0, 9), (0, 10));

        arq.modify(1, 3, &10);
        arq.modify(3, 5, &1);

        assert_eq!(arq.query(0, 9), (23, 10));
        assert_eq!(arq.query(10, 4), (0, 0));
    }

    #[test]
    fn test_dynamic_range_sum() {
        let initializer = Box::new(|l, r| (0, 1 + r - l));
        let mut arq = DynamicArq::<AssignSum>::new(0, 9, false, initializer);

        assert_eq!(arq.query(0, 0, 9), (0, 10));

        arq.modify(0, 1, 3, &10);
        arq.modify(0, 3, 5, &1);

        assert_eq!(arq.query(0, 0, 9), (23, 10));
        assert_eq!(arq.query(0, 10, 4), (0, 0));
    }

    #[test]
    fn test_supply_demand() {
        let mut arq = StaticArq::<SupplyDemand>::new(vec![(0, 0, 0); 10]);

        arq.modify(1, 1, &(25, 100));
        arq.modify(3, 3, &(100, 30));
        arq.modify(9, 9, &(0, 20));

        assert_eq!(arq.query(0, 9), (125, 150, 75));
    }

    #[test]
    fn test_dynamic_supply_demand() {
        let mut arq = DynamicArq::<SupplyDemand>::new_with_identity(0, 9, false);

        arq.modify(0, 1, 1, &(25, 100));
        arq.modify(0, 3, 3, &(100, 30));
        arq.modify(0, 9, 9, &(0, 20));

        assert_eq!(arq.query(0, 0, 9), (125, 150, 75));
    }

    #[test]
    fn test_binary_search_rmq() {
        let vec = vec![2, 1, 0, -1, -2, -3, -4, -5];
        let mut arq = StaticArq::<AssignMin>::new(vec);
        let pos = first_negative_static(&mut arq);

        arq.modify(3, 7, &0);
        let pos_zeros = first_negative_static(&mut arq);

        assert_eq!(pos, 3);
        assert_eq!(pos_zeros, -1);
    }

    #[test]
    fn test_dynamic_binary_search_rmq() {
        let initializer = Box::new(|_, r| 2 - r);
        let (l_bound, r_bound) = (0, 7);
        let mut arq = DynamicArq::<AssignMin>::new(l_bound, r_bound, false, initializer);
        let pos = first_negative_dynamic(&mut arq, 0, l_bound, r_bound);

        arq.modify(0, 2, 7, &0);
        let pos_zeros = first_negative_dynamic(&mut arq, 0, l_bound, r_bound);

        assert_eq!(pos, 3);
        assert_eq!(pos_zeros, -1);
    }

    #[test]
    fn test_mos_algorithm() {
        let queries = vec![(0, 2, ()), (5, 5, ()), (2, 6, ()), (0, 6, ())];
        let arr = vec![4, 8, 4, 7, 1, 9, 8];

        let answers = DistinctVals::new(arr).process(&queries);

        assert_eq!(answers, vec![2, 1, 5, 5]);
    }
}
