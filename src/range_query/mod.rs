pub mod dynamic_arq;
pub mod specs;
pub mod sqrt_decomp;
pub mod static_arq;
pub use dynamic_arq::{ArqView, DynamicArq};
pub use specs::ArqSpec;
pub use static_arq::StaticArq;

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

#[cfg(test)]
mod test {
    use super::specs::*;
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
    fn test_rmq() {
        let mut arq = StaticArq::<AssignMin>::new(&[0; 10]);

        assert_eq!(arq.query(0, 9), 0);

        arq.update(2, 4, &-5);
        arq.update(5, 7, &-3);
        arq.update(1, 6, &1);

        assert_eq!(arq.query(0, 9), -3);
    }

    #[test]
    fn test_dynamic_rmq() {
        let mut arq = DynamicArq::<AssignMin>::new(false);
        let view = arq.build_from_slice(&[0; 10]);

        assert_eq!(arq.query(view, 0, 9), 0);

        arq.update(view, 2, 4, &-5);
        arq.update(view, 5, 7, &-3);
        arq.update(view, 1, 6, &1);

        assert_eq!(arq.query(view, 0, 9), -3);
    }

    #[test]
    fn test_persistent_rmq() {
        let mut arq = DynamicArq::<AssignMin>::new(true);
        let mut view = arq.build_from_slice(&[0; 10]);

        let at_init = view;
        view = arq.update(view, 2, 4, &-5);
        let snapshot = view;
        view = arq.update(view, 5, 7, &-3);
        view = arq.update(view, 1, 6, &1);

        assert_eq!(arq.query(at_init, 0, 9), 0);
        assert_eq!(arq.query(snapshot, 0, 9), -5);
        assert_eq!(arq.query(view, 0, 9), -3);
    }

    #[test]
    fn test_huge_rmq() {
        let quintillion = 1_000_000_000_000_000_000;
        let mut arq = DynamicArq::<AssignMin>::new(false);
        let view = arq.build_from_identity(9 * quintillion + 1);

        arq.update(view, 2 * quintillion, 4 * quintillion, &-5);
        arq.update(view, 5 * quintillion, 7 * quintillion, &-3);
        arq.update(view, 1 * quintillion, 6 * quintillion, &1);

        assert_eq!(arq.query(view, 0, 9 * quintillion), -3);
    }

    #[test]
    fn test_range_sum() {
        let mut arq = StaticArq::<AssignSum>::new(&[0; 10]);

        assert_eq!(arq.query(0, 9), 0);

        arq.update(1, 3, &10);
        arq.update(3, 5, &1);

        assert_eq!(arq.query(0, 9), 23);
        assert_eq!(arq.query(10, 4), 0);
    }

    #[test]
    fn test_dynamic_range_sum() {
        let mut arq = DynamicArq::<AssignSum>::new(false);
        let view = arq.build_from_slice(&[0; 10]);

        assert_eq!(arq.query(view, 0, 9), 0);

        arq.update(view, 1, 3, &10);
        arq.update(view, 3, 5, &1);

        assert_eq!(arq.query(view, 0, 9), 23);
        assert_eq!(arq.query(view, 10, 4), 0);
    }

    #[test]
    fn test_supply_demand() {
        let mut arq = StaticArq::<SupplyDemand>::new(&[(0, 0, 0); 10]);

        arq.update(1, 1, &(25, 100));
        arq.update(3, 3, &(100, 30));
        arq.update(9, 9, &(0, 20));

        assert_eq!(arq.query(0, 9), (125, 150, 75));
    }

    #[test]
    fn test_dynamic_supply_demand() {
        let mut arq = DynamicArq::<SupplyDemand>::new(false);
        let view = arq.build_from_identity(10);

        arq.update(view, 1, 1, &(25, 100));
        arq.update(view, 3, 3, &(100, 30));
        arq.update(view, 9, 9, &(0, 20));

        assert_eq!(arq.query(view, 0, 9), (125, 150, 75));
    }

    #[test]
    fn test_binary_search_rmq() {
        let vec = vec![2, 1, 0, -1, -2, -3, -4, -5];
        let mut arq = StaticArq::<AssignMin>::new(&vec);
        let first_neg = static_arq::first_negative(&mut arq);

        arq.update(3, 7, &0);
        let first_neg_zeros = static_arq::first_negative(&mut arq);

        assert_eq!(first_neg, Some(3));
        assert_eq!(first_neg_zeros, None);
    }

    #[test]
    fn test_dynamic_binary_search_rmq() {
        let vec = vec![2, 1, 0, -1, -2, -3, -4, -5];
        let mut arq = DynamicArq::<AssignMin>::new(false);
        let view = arq.build_from_slice(&vec);
        let first_neg = dynamic_arq::first_negative(&mut arq, view);

        arq.update(view, 3, 7, &0);
        let first_neg_zeros = dynamic_arq::first_negative(&mut arq, view);

        assert_eq!(first_neg, Some(3));
        assert_eq!(first_neg_zeros, None);
    }
}
