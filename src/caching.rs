//! Basic Cacher struct which stores a closure and a hashmap.
//! The hashmap stores key value pairs representing previous
//! function calls.
//!
//! When the Cacher function is run, it first does a lookup
//! to see if the value has already been calculated. If it has,
//! it returns that value. If it hasn't, it calculates the value,
//! adds it to the hashmap, and returns it.

use std::collections::HashMap;

/// The Cacher struct (Memoization) stores a function and a Hashmap.
/// The HashMap keeps track of previous input and output for the function so
/// that it only ever has to be called once per input. Use for expensive functions.
pub struct Cacher<F, U, V>
where
    F: Fn(U) -> V,
    U: std::cmp::Eq + std::hash::Hash + Copy,
    V: Copy,
{
    calculation: F,
    values: HashMap<U, V>,
}

impl<F, U, V> Cacher<F, U, V>
where
    F: Fn(U) -> V,
    U: std::cmp::Eq + std::hash::Hash + Copy,
    V: Copy,
{
    /// Constuctor for the Casher
    /// # Examples
    /// ```
    /// # use contest_algorithms::caching::Cacher;
    /// let mut squared = Cacher::new(|n: u32| n*n);
    /// ```
    pub fn new(calculation: F) -> Cacher<F, U, V> {
        Cacher {
            calculation,
            values: HashMap::new(),
        }
    }

    /// Performs a lookup into the HashMap to see if the value has already
    /// been calculated. If it has, returns the value. If it has not,
    /// calls the function, stores the value, then returns the value.
    /// # Examples
    /// ```
    /// # use contest_algorithms::caching::Cacher;
    /// let mut squared = Cacher::new(|n: u32| n*n);
    ///
    /// // This is where we call the function
    /// let sixteen = squared.call(4);
    /// ```
    // TODO: whenever Rust's Entry API gains the ability to take ownership of
    // arg only when necessary, this method should follow the same practice.
    // Also, Cacher should implement Fn(U)->V once this is possible.
    pub fn call(&mut self, arg: U) -> V {
        let calc = &self.calculation;
        *self.values.entry(arg).or_insert_with_key(|&key| calc(key))
    }

    /// Calls the function without performing a lookup and replaces
    /// the old return value with the new one, and returns it.
    /// Potentially useful if the function reads from a file or RNG
    /// whose state may have changed.
    // TODO: if there's state, FnMut seems more appropriate.
    pub fn call_and_replace(&mut self, arg: U) -> V {
        let new_val = (self.calculation)(arg);
        self.values.insert(arg, new_val);
        new_val
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cacher_basically_works() {
        let mut word_len = Cacher::new(|word: &str| word.len());
        let hello = word_len.call("hello");

        // Test function returns correctly
        assert_eq!(hello, 5);

        // Test HashMap is correct length
        assert_eq!(word_len.values.len(), 1);

        // Test HashMap has correct value after one insert
        let mut test_map = HashMap::new();
        test_map.insert("hello", 5);
        assert_eq!(word_len.values, test_map);

        // Test HashMap has correct value after duplicate insert
        word_len.call("hello");
        assert_eq!(word_len.values, test_map);

        // Test HashMap has correct values after unique input
        word_len.call("wazzup");
        test_map.insert("wazzup", 6);
        assert_eq!(word_len.values, test_map);
    }

    #[test]
    fn test_cacher_speed() {
        // Simulate a function that takes 1 second to complete
        let mut func = Cacher::new(|x| {
            std::thread::sleep(std::time::Duration::from_millis(100));
            x * x
        });

        // Would take 10 minutes without caching
        for _ in 0..6000 {
            assert_eq!(25, func.call(5));
        }
    }

    #[test]
    fn test_call_and_replace() {
        use std::time::Instant;

        let mut func = Cacher::new(|_param: usize| Instant::now());
        let first_instant = func.call(0);
        let lookup_instant = func.call(0);

        assert_eq!(first_instant, lookup_instant);
        assert_eq!(1, func.values.len());

        let second_instant = func.call_and_replace(0);
        assert_eq!(1, func.values.len());
        assert_ne!(second_instant, lookup_instant);
    }
}
