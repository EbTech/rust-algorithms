//! Basic Cacher struct which stores a closure and a hashmap.
//! The hasmap stores key value pairs representing previous
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
    /// calls the function, stores the value, then returns the value
    /// # Examples
    /// ```
    /// # use contest_algorithms::caching::Cacher;
    /// let mut squared = Cacher::new(|n: u32| n*n);
    ///
    /// // This is where we call the function
    /// let sixteen = squared.call(4);
    /// ```
    pub fn call(&mut self, arg: U) -> V {
        // This is basically the magic of the whole
        // structure. You can do this with the entry
        // api, but I like how readable this particular
        // block of code is.
        if let Some(&val) = self.values.get(&arg) {
            val
        } else {
            let val = (self.calculation)(arg);
            self.values.insert(arg, val);
            val
        }
    }

    /// Calls the function without performing a lookup and replaces
    /// the old calculation with the new one, then returns the value
    ///
    /// # Use Case
    /// If you're wondering, this is for if some sort of "state" has changed
    /// underneath you, so your same function call with the same input
    /// might now have different output. For instance, if part of your function
    /// reads from a file and
    /// you think the contents of that file have changed even though the name
    /// has not.
    pub fn call_and_replace(&mut self, arg: U) -> V {
        let new_val = (self.calculation)(arg);
        self.values.insert(arg, new_val);
        new_val
    }
}

#[cfg(test)]
mod tests {

    use super::Cacher;
    use std::collections::HashMap;

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
    fn call_and_replace() {
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
