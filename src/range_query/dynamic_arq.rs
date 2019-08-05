//! Associative Range Query Tree with dynamic allocation, supporting dynamic
//! node construction and persistence
use super::ArqSpec;

pub struct DynamicArqNode<T: ArqSpec> {
    val: T::M,
    app: Option<T::F>,
    down: (usize, usize),
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
            down: (usize::max_value(), usize::max_value()),
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

pub type ArqView = (usize, i64, i64);

/// A dynamic, and optionally persistent, associative range query data structure.
pub struct DynamicArq<T: ArqSpec> {
    nodes: Vec<DynamicArqNode<T>>,
    is_persistent: bool,
    initializer: Box<dyn Fn(i64, i64) -> T::M>,
}

impl<T: ArqSpec> DynamicArq<T> {
    /// Initializes the data structure without creating any nodes.
    /// The initializer f must satisfy: f(l, r) = T::op(f(l, m), f(m+1, r))
    pub fn new(is_persistent: bool, initializer: Box<dyn Fn(i64, i64) -> T::M>) -> Self {
        Self {
            nodes: vec![],
            is_persistent,
            initializer,
        }
    }

    /// Initializes the data structure without creating any nodes.
    pub fn new_with_identity(is_persistent: bool) -> Self {
        let initializer = Box::new(|_, _| T::identity());
        Self::new(is_persistent, initializer)
    }

    /// Lazily builds a tree over the range [l_bound, r_bound].
    pub fn build_using_initializer(&mut self, l_bound: i64, r_bound: i64) -> ArqView {
        let view = (self.nodes.len(), l_bound, r_bound);
        let val = (self.initializer)(l_bound, r_bound);
        self.nodes.push(DynamicArqNode::new(val));
        view
    }

    /// Builds a tree whose leaves are set to a given non-empty slice.
    pub fn build_from_slice(&mut self, init_val: &[T::M]) -> ArqView {
        if init_val.len() == 1 {
            let view = (self.nodes.len(), 0, 0);
            let val = T::op(&T::identity(), &init_val[0]);
            self.nodes.push(DynamicArqNode::new(val));
            return view;
        }
        let m = (init_val.len() + 1) / 2;
        let (l_init, r_init) = init_val.split_at(m);
        let l_view = self.build_from_slice(l_init);
        let r_view = self.build_from_slice(r_init);
        self.merge_equal_sized(l_view, r_view)
    }

    /// Merges two balanced subtrees into a single tree with a 0-indexed view.
    pub fn merge_equal_sized(&mut self, l_view: ArqView, r_view: ArqView) -> ArqView {
        let l_len = l_view.2 - l_view.1 + 1;
        let r_len = r_view.2 - r_view.1 + 1;
        assert!(l_len == r_len || l_len == r_len + 1);

        let view = (self.nodes.len(), 0, l_len + r_len - 1);
        let root = DynamicArqNode {
            val: T::identity(),
            app: None,
            down: (l_view.0, r_view.0),
        };
        self.nodes.push(root);
        self.pull(view.0);
        view
    }

    pub fn push(&mut self, (p, l, r): ArqView) -> (ArqView, ArqView) {
        let m = (l + r) / 2;
        if self.nodes[p].down.0 == usize::max_value() {
            let l_val = (self.initializer)(l, m);
            let r_val = (self.initializer)(m + 1, r);
            self.nodes.push(DynamicArqNode::new(l_val));
            self.nodes.push(DynamicArqNode::new(r_val));
            self.nodes[p].down = (self.nodes.len() - 2, self.nodes.len() - 1)
        };
        let (lp, rp) = self.nodes[p].down;
        if let Some(ref f) = self.nodes[p].app.take() {
            self.nodes[lp].apply(f, l == m);
            self.nodes[rp].apply(f, m + 1 == r);
        }
        ((lp, l, m), (rp, m + 1, r))
    }

    fn pull(&mut self, p: usize) {
        let (lp, rp) = self.nodes[p].down;
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

    /// Applies the endomorphism f to all entries from l to r, inclusive.
    /// If l == r, the updates are eager. Otherwise, they are lazy.
    pub fn modify(&mut self, view: ArqView, l: i64, r: i64, f: &T::F) -> ArqView {
        let (p, cl, cr) = view;
        if r < cl || cr < l {
            view
        } else if l <= cl && cr <= r /* && self.l == self.r forces eager */ {
            let p_clone = self.clone_node(p);
            self.nodes[p_clone].apply(f, l == r);
            (p_clone, cl, cr)
        } else {
            let (l_view, r_view) = self.push(view);
            let p_clone = self.clone_node(p);
            let lp_clone = self.modify(l_view, l, r, f).0;
            let rp_clone = self.modify(r_view, l, r, f).0;
            self.nodes[p_clone].down = (lp_clone, rp_clone);
            self.pull(p_clone);
            (p_clone, cl, cr)
        }
    }

    /// Returns the aggregate range query on all entries from l to r, inclusive.
    pub fn query(&mut self, view: ArqView, l: i64, r: i64) -> T::M {
        let (p, cl, cr) = view;
        if r < cl || cr < l {
            T::identity()
        } else if l <= cl && cr <= r {
            T::op(&T::identity(), &self.nodes[p].val)
        } else {
            let (l_view, r_view) = self.push(view);
            let l_agg = self.query(l_view, l, r);
            let r_agg = self.query(r_view, l, r);
            T::op(&l_agg, &r_agg)
        }
    }
}

/// An example of binary search on the tree of a DynamicArq.
/// The tree may have any size, not necessarily a power of two.
pub fn first_negative(arq: &mut DynamicArq<super::specs::AssignMin>, view: ArqView) -> Option<i64> {
    let (p, cl, cr) = view;
    if cl == cr {
        Some(cl).filter(|_| arq.nodes[p].val < 0)
    } else {
        let (l_view, r_view) = arq.push(view);
        let lp = l_view.0;
        if arq.nodes[lp].val < 0 {
            first_negative(arq, l_view)
        } else {
            first_negative(arq, r_view)
        }
    }
}
