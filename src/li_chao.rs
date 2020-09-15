/// A structure for answering minimum queries on a set of linear functions. Supports two
/// operations: inserting a linear function and querying for minimum at a given point. Unlike the
/// simplest convex hull trick implementation, the queries can be done in any order, and we can do
/// all the calculations using integers.

/// This data structure builds a static segment tree over the interval of x-coordinates
/// [left,right), so requires O(right-left) memory. Just like normal segment trees, this could be
/// modified to a dynamic tree when the range is huge (it can also be made persistent!).
/// It can also probably be done as an AlCash style iterative segment tree in the static case.

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
            lines: vec![(0, i64::MAX); 4 * (right - left + 1) as usize],
        }
    }

    /// Every node in the tree has the property that the line that minimizes its midpoint is found
    /// either in the node or one of its ancestors.  When we visit a node, we compute the winner at
    /// the midpoint of the node. The winner is stored in the node. The loser can still possibly
    /// beat the winner on some segment, either to the left or to the right of the current
    /// midpoint, so we propagate it to that segment. This sequence ensures that the invariant is
    /// kept.
    fn add_line_impl(&mut self, mut m: i64, mut b: i64, ix: usize, l: i64, r: i64) {
        let x = (l + r) / 2;
        if m * x + b < self.lines[ix].0 * x + self.lines[ix].1 {
            std::mem::swap(&mut m, &mut self.lines[ix].0);
            std::mem::swap(&mut b, &mut self.lines[ix].1);
        }
        if r - l > 1 {
            if m < self.lines[ix].0 {
                self.add_line_impl(m, b, 2 * ix + 1, x, r);
            } else {
                self.add_line_impl(m, b, 2 * ix, l, x);
            }
        }
    }

    /// Adds the line with slope m and intercept b. O(log N) complexity.
    pub fn add_line(&mut self, m: i64, b: i64) {
        self.add_line_impl(m, b, 1, self.left, self.right);
    }

    /// Because of the invariant established by add_line, we know that the best line for a given
    /// point is stored in one of the ancestors of its node. So we accumulate the minimum answer as
    /// we go back up the tree.
    fn query_impl(&self, x: i64, ix: usize, l: i64, r: i64) -> i64 {
        let y = self.lines[ix].0 * x + self.lines[ix].1;
        if r - l == 1 {
            y
        } else if x >= (l + r) / 2 {
            self.query_impl(x, 2 * ix + 1, (l + r) / 2, r).min(y)
        } else {
            self.query_impl(x, 2 * ix, l, (l + r) / 2).min(y)
        }
    }

    /// Finds the minimum mx+b among all lines in the structure. O(log N) complexity.
    pub fn query(&self, x: i64) -> i64 {
        self.query_impl(x, 1, self.left, self.right)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_li_chao_tree() {
        let lines = [(0, 3), (1, 0), (-1, 8), (2, -1), (-1, 4)];
        let xs = [0, 1, 2, 3, 4, 5];
        // results[i] consists of the expected y-coordinates after processing
        // the first i+1 lines.
        let results = [
            [3, 3, 3, 3, 3, 3],
            [0, 1, 2, 3, 3, 3],
            [0, 1, 2, 3, 3, 3],
            [-1, 1, 2, 3, 3, 3],
            [-1, 1, 2, 1, 0, -1],
        ];
        let mut li_chao = LiChaoTree::new(0, 6);

        assert_eq!(li_chao.query(0), i64::MAX);
        for (&(slope, intercept), expected) in lines.iter().zip(results.iter()) {
            li_chao.add_line(slope, intercept);
            let ys: Vec<i64> = xs.iter().map(|&x| li_chao.query(x)).collect();
            assert_eq!(expected, &ys[..]);
        }
    }
}
