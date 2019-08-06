pub mod dynamic_arq;
pub mod specs;
pub mod sqrt_decomp;
pub mod static_arq;
pub use dynamic_arq::{ArqView, DynamicArq};
pub use specs::ArqSpec;
pub use static_arq::StaticArq;

#[cfg(test)]
mod test {
    use super::specs::*;
    use super::*;

    #[test]
    fn test_rmq() {
        let mut arq = StaticArq::<AssignMin>::new(&[0; 10]);

        assert_eq!(arq.query(0, 9), 0);

        arq.modify(2, 4, &-5);
        arq.modify(5, 7, &-3);
        arq.modify(1, 6, &1);

        assert_eq!(arq.query(0, 9), -3);
    }

    #[test]
    fn test_dynamic_rmq() {
        let mut arq = DynamicArq::<AssignMin>::new(false);
        let view = arq.build_from_slice(&[0; 10]);

        assert_eq!(arq.query(view, 0, 9), 0);

        arq.modify(view, 2, 4, &-5);
        arq.modify(view, 5, 7, &-3);
        arq.modify(view, 1, 6, &1);

        assert_eq!(arq.query(view, 0, 9), -3);
    }

    #[test]
    fn test_persistent_rmq() {
        let mut arq = DynamicArq::<AssignMin>::new(true);
        let mut view = arq.build_from_slice(&[0; 10]);

        let at_init = view;
        view = arq.modify(view, 2, 4, &-5);
        let snapshot = view;
        view = arq.modify(view, 5, 7, &-3);
        view = arq.modify(view, 1, 6, &1);

        assert_eq!(arq.query(at_init, 0, 9), 0);
        assert_eq!(arq.query(snapshot, 0, 9), -5);
        assert_eq!(arq.query(view, 0, 9), -3);
    }

    #[test]
    fn test_huge_rmq() {
        let quintillion = 1_000_000_000_000_000_000;
        let mut arq = DynamicArq::<AssignMin>::new(false);
        let view = arq.build_from_identity(9 * quintillion + 1);

        arq.modify(view, 2 * quintillion, 4 * quintillion, &-5);
        arq.modify(view, 5 * quintillion, 7 * quintillion, &-3);
        arq.modify(view, 1 * quintillion, 6 * quintillion, &1);

        assert_eq!(arq.query(view, 0, 9 * quintillion), -3);
    }

    #[test]
    fn test_range_sum() {
        let mut arq = StaticArq::<AssignSum>::new(&[(0, 1); 10]);

        assert_eq!(arq.query(0, 9), (0, 10));

        arq.modify(1, 3, &10);
        arq.modify(3, 5, &1);

        assert_eq!(arq.query(0, 9), (23, 10));
        assert_eq!(arq.query(10, 4), (0, 0));
    }

    #[test]
    fn test_dynamic_range_sum() {
        let mut arq = DynamicArq::<AssignSum>::new(false);
        let view = arq.build_from_slice(&[(0, 1); 10]);

        assert_eq!(arq.query(view, 0, 9), (0, 10));

        arq.modify(view, 1, 3, &10);
        arq.modify(view, 3, 5, &1);

        assert_eq!(arq.query(view, 0, 9), (23, 10));
        assert_eq!(arq.query(view, 10, 4), (0, 0));
    }

    #[test]
    fn test_supply_demand() {
        let mut arq = StaticArq::<SupplyDemand>::new(&[(0, 0, 0); 10]);

        arq.modify(1, 1, &(25, 100));
        arq.modify(3, 3, &(100, 30));
        arq.modify(9, 9, &(0, 20));

        assert_eq!(arq.query(0, 9), (125, 150, 75));
    }

    #[test]
    fn test_dynamic_supply_demand() {
        let mut arq = DynamicArq::<SupplyDemand>::new(false);
        let view = arq.build_from_identity(10);

        arq.modify(view, 1, 1, &(25, 100));
        arq.modify(view, 3, 3, &(100, 30));
        arq.modify(view, 9, 9, &(0, 20));

        assert_eq!(arq.query(view, 0, 9), (125, 150, 75));
    }

    #[test]
    fn test_binary_search_rmq() {
        let vec = vec![2, 1, 0, -1, -2, -3, -4, -5];
        let mut arq = StaticArq::<AssignMin>::new(&vec);
        let first_neg = static_arq::first_negative(&mut arq);

        arq.modify(3, 7, &0);
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

        arq.modify(view, 3, 7, &0);
        let first_neg_zeros = dynamic_arq::first_negative(&mut arq, view);

        assert_eq!(first_neg, Some(3));
        assert_eq!(first_neg_zeros, None);
    }
}
