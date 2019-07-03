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

impl<T: ArqSpec> ArqTree<T>
where
    T::F: Clone,
{
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
        for s in (1..32).rev() {
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
    /// Panics if l or r is out of range.
    pub fn modify(&mut self, mut l: usize, mut r: usize, f: &T::F) {
        l += self.app.len();
        r += self.app.len();
        self.push_to(l);
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
    /// Panics if l or r is out of range.
    pub fn query(&mut self, mut l: usize, mut r: usize) -> T::M {
        l += self.app.len();
        r += self.app.len();
        self.push_to(l);
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
    type F;
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
pub struct AssignMin;
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
pub struct AssignSum;
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
pub struct SupplyDemand;
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
    }

    #[test]
    fn test_supply_demand() {
        let mut arq = ArqTree::<SupplyDemand>::new(vec![(0, 0, 0); 10]);

        arq.modify(1, 1, &(25, 100));
        arq.modify(3, 3, &(100, 30));
        arq.modify(9, 9, &(0, 20));

        assert_eq!(arq.query(0, 9), (125, 150, 75));
    }
}
