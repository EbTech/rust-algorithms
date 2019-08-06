//! A collection of example ArqSpec implementations

pub trait ArqSpec {
    /// Type of data representing an endomorphism.
    // Note that while a Fn(M) -> M may seem like a more natural representation
    // for an endomorphism, compositions would then have to delegate to each of
    // their parts. This representation is more efficient.
    type F: Clone;
    /// Type of monoid elements.
    type M: Clone;

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
        i64::max_value()
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
