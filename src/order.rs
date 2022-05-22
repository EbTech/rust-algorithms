//! Ordering algorithms.

/// A comparator on partially ordered elements, that panics if they are incomparable
///
/// # Example
///
/// ```
/// use contest_algorithms::order::asserting_cmp;
/// let mut vec = vec![4.5, -1.7, 1.2];
/// vec.sort_unstable_by(asserting_cmp);
/// assert_eq!(vec, vec![-1.7, 1.2, 4.5]);
/// ```
pub fn asserting_cmp<T: PartialOrd>(a: &T, b: &T) -> std::cmp::Ordering {
    a.partial_cmp(b).expect("Comparing incomparable elements")
}

/// Assuming slice is sorted and totally ordered, returns the minimum i for which
/// slice[i] >= key, or slice.len() if no such i exists
pub fn slice_lower_bound<T: PartialOrd>(slice: &[T], key: &T) -> usize {
    slice
        .binary_search_by(|x| asserting_cmp(x, key).then(std::cmp::Ordering::Greater))
        .unwrap_err()
}

/// Assuming slice is sorted and totally ordered, returns the minimum i for which
/// slice[i] > key, or slice.len() if no such i exists
pub fn slice_upper_bound<T: PartialOrd>(slice: &[T], key: &T) -> usize {
    slice
        .binary_search_by(|x| asserting_cmp(x, key).then(std::cmp::Ordering::Less))
        .unwrap_err()
}

/// Stably merges two sorted and totally ordered collections into one
pub fn merge_sorted<T: PartialOrd>(
    i1: impl IntoIterator<Item = T>,
    i2: impl IntoIterator<Item = T>,
) -> Vec<T> {
    let mut i1 = i1.into_iter().peekable();
    let mut i2 = i2.into_iter().peekable();
    let mut merged = Vec::with_capacity(i1.size_hint().0 + i2.size_hint().0);
    while let (Some(a), Some(b)) = (i1.peek(), i2.peek()) {
        merged.push(if a <= b { i1.next() } else { i2.next() }.unwrap());
    }
    merged.extend(i1.chain(i2));
    merged
}

/// A stable sort
pub fn merge_sort<T: Ord>(mut v: Vec<T>) -> Vec<T> {
    if v.len() < 2 {
        v
    } else {
        let v2 = v.split_off(v.len() / 2);
        merge_sorted(merge_sort(v), merge_sort(v2))
    }
}

/// A simple data structure for coordinate compression
pub struct SparseIndex {
    coords: Vec<i64>,
}

impl SparseIndex {
    /// Builds an index, given the full set of coordinates to compress.
    pub fn new(mut coords: Vec<i64>) -> Self {
        coords.sort_unstable();
        coords.dedup();
        Self { coords }
    }

    /// Returns Ok(i) if the coordinate q appears at index i
    /// Returns Err(i) if q appears between indices i-1 and i
    pub fn compress(&self, q: i64) -> Result<usize, usize> {
        self.coords.binary_search(&q)
    }
}

/// Represents a maximum (upper envelope) of a collection of linear functions of one
/// variable, evaluated using an online version of the convex hull trick.
/// It combines the offline algorithm with square root decomposition, resulting in an
/// asymptotically suboptimal but simple algorithm with good amortized performance:
/// N inserts interleaved with Q queries yields O(N sqrt Q + Q log N) time complexity
/// in general, or O((N + Q) log N) if all queries come after all inserts.
// Proof: the Q log N term comes from calls to slice_lower_bound(). As for the N sqrt Q,
//        note that between successive times when the hull is rebuilt, O(N) work is done,
//        and the running totals of insertions and queries satisfy del_N (del_Q + 1) > N.
//        Now, either del_Q >= sqrt Q, or else del_Q <= 2 sqrt Q - 1
//                                          => del_N > N / (2 sqrt Q).
//        Since del(N sqrt Q) >= max(N del(sqrt Q), del_N sqrt Q)
//                            >= max(N del_Q / (2 sqrt Q), del_N sqrt Q),
//        we conclude that del(N sqrt Q) >= N / 2.
#[derive(Default)]
pub struct PiecewiseLinearConvexFn {
    recent_lines: Vec<(f64, f64)>,
    sorted_lines: Vec<(f64, f64)>,
    intersections: Vec<f64>,
    amortized_work: usize,
}

impl PiecewiseLinearConvexFn {
    /// Replaces the represented function with the maximum of itself and a provided line
    pub fn max_with(&mut self, new_m: f64, new_b: f64) {
        self.recent_lines.push((new_m, new_b));
    }

    /// Similar to max_with but requires that (new_m, new_b) be the largest pair so far
    fn max_with_sorted(&mut self, new_m: f64, new_b: f64) {
        while let Some(&(last_m, last_b)) = self.sorted_lines.last() {
            // If slopes are equal, get rid of the old line as its intercept is lower
            if (new_m - last_m).abs() > 1e-9 {
                let intersect = (new_b - last_b) / (last_m - new_m);
                if self.intersections.last() < Some(&intersect) {
                    self.intersections.push(intersect);
                    break;
                }
            }
            self.intersections.pop();
            self.sorted_lines.pop();
        }
        self.sorted_lines.push((new_m, new_b));
    }

    /// Evaluates the function at x
    fn eval_unoptimized(&self, x: f64) -> f64 {
        let idx = slice_lower_bound(&self.intersections, &x);
        self.recent_lines
            .iter()
            .chain(self.sorted_lines.get(idx))
            .map(|&(m, b)| m * x + b)
            .max_by(asserting_cmp)
            .unwrap_or(-1e18)
    }

    /// Evaluates the function at x with good amortized runtime
    pub fn evaluate(&mut self, x: f64) -> f64 {
        self.amortized_work += self.recent_lines.len();
        if self.amortized_work > self.sorted_lines.len() {
            self.amortized_work = 0;
            self.recent_lines.sort_unstable_by(asserting_cmp);
            self.intersections.clear();
            let all_lines = merge_sorted(self.recent_lines.drain(..), self.sorted_lines.drain(..));
            for (new_m, new_b) in all_lines {
                self.max_with_sorted(new_m, new_b);
            }
        }
        self.eval_unoptimized(x)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bounds() {
        let mut vals = vec![16, 45, 45, 45, 82];

        assert_eq!(slice_upper_bound(&vals, &44), 1);
        assert_eq!(slice_lower_bound(&vals, &45), 1);
        assert_eq!(slice_upper_bound(&vals, &45), 4);
        assert_eq!(slice_lower_bound(&vals, &46), 4);

        vals.dedup();
        for (i, q) in vals.iter().enumerate() {
            assert_eq!(slice_lower_bound(&vals, q), i);
            assert_eq!(slice_upper_bound(&vals, q), i + 1);
        }
    }

    #[test]
    fn test_merge_sorted() {
        let vals1 = vec![16, 45, 45, 82];
        let vals2 = vec![-20, 40, 45, 50];
        let vals_merged = vec![-20, 16, 40, 45, 45, 45, 50, 82];

        assert_eq!(merge_sorted(None, Some(42)), vec![42]);
        assert_eq!(merge_sorted(vals1.iter().cloned(), None), vals1);
        assert_eq!(merge_sorted(vals1, vals2), vals_merged);
    }

    #[test]
    fn test_merge_sort() {
        let unsorted = vec![8, -5, 1, 4, -3, 4];
        let sorted = vec![-5, -3, 1, 4, 4, 8];

        assert_eq!(merge_sort(unsorted), sorted);
        assert_eq!(merge_sort(sorted.clone()), sorted);
    }

    #[test]
    fn test_coord_compress() {
        let mut coords = vec![16, 99, 45, 18];
        let index = SparseIndex::new(coords.clone());

        coords.sort_unstable();
        for (i, q) in coords.into_iter().enumerate() {
            assert_eq!(index.compress(q - 1), Err(i));
            assert_eq!(index.compress(q), Ok(i));
            assert_eq!(index.compress(q + 1), Err(i + 1));
        }
    }

    #[test]
    fn test_range_compress() {
        let queries = vec![(0, 10), (10, 19), (20, 29)];
        let coords = queries.iter().flat_map(|&(i, j)| vec![i, j + 1]).collect();
        let index = SparseIndex::new(coords);

        assert_eq!(index.coords, vec![0, 10, 11, 20, 30]);
    }

    #[test]
    fn test_convex_hull_trick() {
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
        let mut func = PiecewiseLinearConvexFn::default();
        assert_eq!(func.evaluate(0.0), -1e18);
        for (&(slope, intercept), expected) in lines.iter().zip(results.iter()) {
            func.max_with(slope as f64, intercept as f64);
            let ys: Vec<i64> = xs.iter().map(|&x| func.evaluate(x as f64) as i64).collect();
            assert_eq!(expected, &ys[..]);
        }
    }
}
