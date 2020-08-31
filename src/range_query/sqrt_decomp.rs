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

/// Represents a minimum (lower envelope) of a collection of linear functions of a variable,
/// evaluated using the convex hull trick with square root decomposition.
#[derive(Debug)]
pub struct PiecewiseLinearFn {
    sorted_lines: Vec<(i64, i64)>,
    intersections: Vec<f64>,
    recent_lines: Vec<(i64, i64)>,
    merge_threshold: usize,
}

impl PiecewiseLinearFn {
    /// For N inserts interleaved with Q queries, a threshold of N/sqrt(Q) yields
    /// O(N sqrt Q + Q log N) time complexity. If all queries come after all inserts,
    /// a threshold of 0 yields O(N + Q log N) time complexity.
    pub fn with_merge_threshold(merge_threshold: usize) -> Self {
        Self {
            sorted_lines: vec![],
            intersections: vec![],
            recent_lines: vec![],
            merge_threshold,
        }
    }

    /// Replaces this function with the minimum of itself and a provided line
    pub fn min_with(&mut self, slope: i64, intercept: i64) {
        self.recent_lines.push((slope, intercept));
    }

    fn update_envelope(&mut self) {
        self.recent_lines.extend(self.sorted_lines.drain(..));
        self.recent_lines.sort_unstable_by_key(|&(m,b)| (-m, b));
        self.intersections.clear();

        'outer: for (m1, b1) in self.recent_lines.drain(..) {
            while let Some(&(m2, b2)) = self.sorted_lines.last() {
                // The lines are sorted by intercept so an incoming line with
                // the same slope can never beat the previous line.
                if m2 == m1 {
                    continue 'outer;
                }
                let new_intersection = (b1 - b2) as f64 / (m2 - m1) as f64;
                if &new_intersection > self.intersections.last().unwrap_or(&f64::MIN) {
                    self.intersections.push(new_intersection);
                    break;
                }
                self.intersections.pop();
                self.sorted_lines.pop();
            }
            self.sorted_lines.push((m1, b1));
        }
    }

    fn eval_in_envelope(&self, x: i64) -> i64 {
        if self.sorted_lines.is_empty() {
            return i64::MAX;
        }
        // Wow is this messy
        let idx = match self
            .intersections
            .binary_search_by(|&y| y.partial_cmp(&(x as f64)).unwrap())
        {
            Ok(k) => k,
            Err(k) => k,
        };
        let (m, b) = self.sorted_lines[idx];
        m * x + b
    }

    fn eval_helper(&self, x: i64) -> i64 {
        self.recent_lines
            .iter()
            .map(|&(m, b)| m * x + b)
            .chain(std::iter::once(self.eval_in_envelope(x)))
            .min()
            .unwrap()
    }

    /// Evaluates the function at x
    pub fn evaluate(&mut self, x: i64) -> i64 {
        if self.recent_lines.len() > self.merge_threshold {
            self.update_envelope();
        }
        self.eval_helper(x)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_mos_algorithm() {
        let queries = vec![(0, 2, ()), (5, 5, ()), (2, 6, ()), (0, 6, ())];
        let arr = vec![4, 8, 4, 7, 1, 9, 8];

        let answers = DistinctVals::new(arr).process(&queries);

        assert_eq!(answers, vec![2, 1, 5, 5]);
    }

    #[test]
    fn test_convex_hull_trick() {
        let lines = [(0, 4), (1, 0), (-1, 10), (2, -1), (-1, 6)];
        let xs = [0, 1, 2, 3, 4, 5];
        let results = [[4, 4, 4, 4, 4, 4],
        [0,1,2,3,4,4]];

        for threshold in 0..=lines.len() {
            let mut func = PiecewiseLinearFn::with_merge_threshold(threshold);
            for (&(slope, intercept), expected) in lines.iter().zip(results.iter()) {
                func.min_with(slope, intercept);
                let ys: Vec<i64> = xs.iter().map(|&x| func.evaluate(x)).collect();
                dbg!(&func, &expected, &ys);
                assert_eq!(expected, &ys[..]);
            }
        }
    }
}
