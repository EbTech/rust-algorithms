//! Basic casher struct which store a closure and a hashmap.
//! The hasmap stores key, value pairs representing previous
//! function calls, so that each time the function is called on an
//! already called value, it simply does a lookup rather than performing
//! a calculation.

mod generic_casher {
    use std::collections::HashMap;

    pub struct Casher<T, U, V>
    where
        T: Fn(U) -> V,
        U: std::cmp::Eq + std::hash::Hash + Copy,
        V: Copy,
    {
        calculation: T,
        values: HashMap<U, V>,
    }

    impl<T, U, V> Casher<T, U, V>
    where
        T: Fn(U) -> V,
        U: std::cmp::Eq + std::hash::Hash + Copy,
        V: Copy,
    {
        pub fn new(calculation: T) -> Casher<T, U, V> {
            Casher {
                calculation,
                values: HashMap::new(),
            }
        }

        pub fn value(&mut self, arg: U) -> V {
            if let Some(val) = self.values.get(&arg) {
                *val
            } else {
                let val = (self.calculation)(arg);
                self.values.insert(arg, val);
                val
            }
        }
    }
}

mod basic_casher_for_readability {
    use std::collections::HashMap;

    pub struct Casher<T> where T: Fn(u32) -> u32 {
        calculation: T,
        values: HashMap<u32,u32>
    }

    impl<T> Casher<T> where T: Fn(u32) -> u32 {
        pub fn new(calculation: T) -> Self {
            Casher {
                calculation,
                values: HashMap::new()
            }
        }

        pub fn value(&mut self, arg: u32) -> u32 {
            if let Some(val) = self.values.get(&arg) {
                *val
            } else {
                let val = (self.calculation)(arg);
                self.values.insert(arg, val);
                val
            }
        }
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn example_usage_basic_casher() {
        use super::basic_casher_for_readability::Casher;

        // Simulating a function that takes 5 seconds to complete
        let mut func_tracker = Casher::new(|x| {
            std::thread::sleep(std::time::Duration::from_secs(5));
            x * x
        });

        // The combination of these tests only takes 5 seconds, not 15
        assert_eq!(25, func_tracker.value(5)); // takes 5 seconds
        assert_eq!(25, func_tracker.value(5)); // gets from hashmap
        assert_eq!(25, func_tracker.value(5)); // gets from hashmap
    }

    #[test]
    fn example_usage_generic_casher() {
        use super::generic_casher::Casher;

        // Simulating a function that takes 5 seconds to complete
        let mut func_tracker = Casher::new(|x: &str| {
            std::thread::sleep(std::time::Duration::from_secs(5));
            x.len()
        });

        // Finishes in 5 seconds rather than 15.
        assert_eq!(5, func_tracker.value("hello"));
        assert_eq!(5, func_tracker.value("hello"));
        assert_eq!(5, func_tracker.value("hello"));
    }
}