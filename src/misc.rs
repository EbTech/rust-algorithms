//! Miscellaneous algorithms.

/// A comparator on partially ordered elements, that panics if they are incomparable
pub fn asserting_cmp<T: PartialOrd>(a: &T, b: &T) -> std::cmp::Ordering {
    a.partial_cmp(b).expect("Comparing incomparable elements")
}

/// Assuming slice is totally ordered and sorted, returns the minimum i for which
/// slice[i] >= key, or slice.len() if no such i exists
pub fn slice_lower_bound<T: PartialOrd>(slice: &[T], key: &T) -> usize {
    slice
        .binary_search_by(|x| asserting_cmp(x, key).then(std::cmp::Ordering::Greater))
        .unwrap_err()
}

/// Assuming slice is totally ordered and sorted, returns the minimum i for which
/// slice[i] > key, or slice.len() if no such i exists
pub fn slice_upper_bound<T: PartialOrd>(slice: &[T], key: &T) -> usize {
    slice
        .binary_search_by(|x| asserting_cmp(x, key).then(std::cmp::Ordering::Less))
        .unwrap_err()
}

/// Merge two sorted collections into one
pub fn merge_sorted<T: PartialOrd>(
    i1: impl IntoIterator<Item = T>,
    i2: impl IntoIterator<Item = T>,
) -> Vec<T> {
    let (mut i1, mut i2) = (i1.into_iter().peekable(), i2.into_iter().peekable());
    let mut merged = Vec::with_capacity(i1.size_hint().0 + i2.size_hint().0);
    while let (Some(a), Some(b)) = (i1.peek(), i2.peek()) {
        merged.push(if a <= b { i1.next() } else { i2.next() }.unwrap());
    }
    merged.extend(i1.chain(i2));
    merged
}

/// A simple data structure for coordinate compression
pub struct SparseIndex {
    coords: Vec<i64>,
}

impl SparseIndex {
    /// Build an index, given the full set of coordinates to compress.
    pub fn new(mut coords: Vec<i64>) -> Self {
        coords.sort_unstable();
        coords.dedup();
        Self { coords }
    }

    /// Return Ok(i) if the coordinate q appears at index i
    /// Return Err(i) if q appears between indices i-1 and i
    pub fn compress(&self, q: i64) -> Result<usize, usize> {
        self.coords.binary_search(&q)
    }
}

/// Represents a maximum (upper envelope) of a collection of linear functions of a variable,
/// evaluated using the convex hull trick with square root decomposition.
pub struct PiecewiseLinearFn {
    sorted_lines: Vec<(f64, f64)>,
    intersections: Vec<f64>,
    recent_lines: Vec<(f64, f64)>,
    merge_threshold: usize,
}

impl PiecewiseLinearFn {
    /// For N inserts interleaved with Q queries, a threshold of N/sqrt(Q) yields
    /// O(N sqrt Q + Q log N) time complexity. If all queries come after all inserts,
    /// any threshold less than N (e.g., 0) yields O(N + Q log N) time complexity.
    pub fn with_merge_threshold(merge_threshold: usize) -> Self {
        Self {
            sorted_lines: vec![],
            intersections: vec![],
            recent_lines: vec![],
            merge_threshold,
        }
    }

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
        if self.recent_lines.len() > self.merge_threshold {
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
        for threshold in 0..=lines.len() {
            let mut func = PiecewiseLinearFn::with_merge_threshold(threshold);
            assert_eq!(func.evaluate(0.0), -1e18);
            for (&(slope, intercept), expected) in lines.iter().zip(results.iter()) {
                func.max_with(slope as f64, intercept as f64);
                let ys: Vec<i64> = xs.iter().map(|&x| func.evaluate(x as f64) as i64).collect();
                assert_eq!(expected, &ys[..]);
            }
        }
    }
}
