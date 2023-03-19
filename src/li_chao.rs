/// A structure for answering maximum queries on a set of linear functions. Supports two
/// operations: inserting a linear function and querying for maximum at a given point. 
/// The queries can be done in any order, and we can do all the calculations using integers.
/// https://cp-algorithms.com/geometry/convex_hull_trick.html#li-chao-tree
/// Compared to the code in the above link, this implementation further improves the algorithm by
/// reducing the number of nodes to (right - left). This is done by removing the midpoint of a
/// segment from both children. Even better, this allows the index of a node to just be the
/// midpoint of the interval!

/// Just like normal segment trees, this could be modified to a dynamic tree when the range is
/// huge, or if the queries are known in advance the x-coordinates can be compressed.
/// (it can also be made persistent!).

pub struct LiChaoTree {
    left: i64,
    right: i64,
    lines: Vec<(i64, i64)>,
}

impl LiChaoTree {
    /// Creates a new tree, built to handle queries on the interval [left, right).
    pub fn new(left: i64, right: i64) -> Self {
        Self {
            left,
            right,
            lines: vec![(0, std::i64::MIN); (right - left) as usize],
        }
    }

    /// Every node in the tree has the property that the line that maximizes its midpoint is found
    /// either in the node or one of its ancestors.  When we visit a node, we compute the winner at
    /// the midpoint of the node. The winner is stored in the node. The loser can still possibly
    /// beat the winner on some segment, either to the left or to the right of the current
    /// midpoint, so we propagate it to that segment. This sequence ensures that the invariant is
    /// kept.
    fn max_with_impl(&mut self, mut m: i64, mut b: i64, l: i64, r: i64) {
        if r <= l {
            return;
        }
        let ix = ((r - self.left + l - self.left) / 2) as usize;
        let mid = self.left + (ix as i64);
        let (ref mut m_ix, ref mut b_ix) = self.lines[ix];
        if m * mid + b > *m_ix * mid + *b_ix {
            std::mem::swap(&mut m, m_ix);
            std::mem::swap(&mut b, b_ix);
        }
        if m < *m_ix {
            self.max_with_impl(m, b, l, mid);
        } else if m > *m_ix {
            self.max_with_impl(m, b, mid + 1, r);
        }
    }

    /// Adds the line with slope m and intercept b. O(log N) complexity.
    pub fn max_with(&mut self, m: i64, b: i64) {
        self.max_with_impl(m, b, self.left, self.right);
    }

    /// Because of the invariant established by add_line, we know that the best line for a given
    /// point is stored in one of the ancestors of its node. So we accumulate the maximum answer as
    /// we go back up the tree.
    fn evaluate_impl(&self, x: i64, l: i64, r: i64) -> i64 {
        if r == l {
            return i64::MIN;
        }
        let ix = ((r - self.left + l - self.left) / 2) as usize;
        let mid = ix as i64 + self.left;
        let y = self.lines[ix].0 * x + self.lines[ix].1;
        if x == mid {
            y
        } else if x < mid {
            self.evaluate_impl(x, l, mid).max(y)
        } else {
            self.evaluate_impl(x, mid + 1, r).max(y)
        }
    }

    /// Finds the maximum mx+b among all lines in the structure. O(log N) complexity.
    pub fn evaluate(&self, x: i64) -> i64 {
        self.evaluate_impl(x, self.left, self.right)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_li_chao_tree() {
        let lines = [(0, -3), (-1, 0), (1, -8), (-2, 1), (1, -4)];
        let xs = [0, 1, 2, 3, 4, 5];
        // results[i] consists of the expected y-coordinates after processing
        // the first i+1 lines.
        let results = [
            [-3, -3, -3, -3, -3, -3],
            [0, -1, -2, -3, -3, -3],
            [0, -1, -2, -3, -3, -3],
            [1, -1, -2, -3, -3, -3],
            [1, -1, -2, -1, 0, 1],
        ];
        let mut li_chao = LiChaoTree::new(0, 6);

        assert_eq!(li_chao.evaluate(0), std::i64::MIN);
        for (&(slope, intercept), expected) in lines.iter().zip(results.iter()) {
            li_chao.max_with(slope, intercept);
            let ys: Vec<i64> = xs.iter().map(|&x| li_chao.evaluate(x)).collect();
            assert_eq!(expected, &ys[..]);
        }
    }
}
