//! Basic casher struct which stores a closure and a hashmap.
//! The hasmap stores key value pairs representing previous
//! function calls.
//!
//! When the Casher function is run, it first does a lookup
//! to see if the value has already been calculated. If it has,
//! it returns that value. If it hasn't, it calculates the value,
//! adds it to the hashmap, and returns it.

use std::collections::HashMap;

/// The Casher struct (Memoization) stores a function and a Hashmap.
/// The HashMap keeps track of previous input and output for the function so
/// that it only ever has to be called once. Use for expensive functions.
pub struct Casher<F, U, V>
where
    F: Fn(U) -> V,
    U: std::cmp::Eq + std::hash::Hash + Copy,
    V: Copy,
{
    calculation: F,
    values: HashMap<U, V>,
}

impl<F, U, V> Casher<F, U, V>
where
    F: Fn(U) -> V,
    U: std::cmp::Eq + std::hash::Hash + Copy,
    V: Copy,
{
    /// Constuctor for the Casher
    /// # Examples
    /// ```
    /// # use contest_algorithms::cashing::Casher;
    /// let mut squared = Casher::new(|n: u32| n*n);
    /// ```
    pub fn new(calculation: F) -> Casher<F, U, V> {
        Casher {
            calculation,
            values: HashMap::new(),
        }
    }

    /// Performs a lookup into the hashmap to see if the value has already
    /// been calculated. If it has, returns the value. If it has not,
    /// calls the function, stores the value, then returns the value
    /// # Examples
    /// ```
    /// # use contest_algorithms::cashing::Casher;
    /// let mut squared = Casher::new(|n: u32| n*n);
    ///
    /// // This is where we call the function
    /// let sixteen = squared.on(4);
    /// ```
    pub fn on(&mut self, arg: U) -> V {
        // Clippy doesn't like that this isn't lazy but I'm fine with it.
        *self.values.entry(arg).or_insert((self.calculation)(arg))
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

    /// Empties the HashMap so that all previous function calls
    /// are forgotten.
    ///
    /// # Examples
    /// ```
    /// # use contest_algorithms::cashing::Casher;
    /// let mut squared = Casher::new(|n: u32| n*n);
    /// let nine = squared.on(3);
    /// let four = squared.on(2);
    ///
    /// // Forget these values
    /// squared.forget();
    /// ```
    pub fn forget(&mut self) {
        self.values = HashMap::new();
    }
}

#[cfg(test)]
mod tests {

    use super::Casher;
    use std::collections::HashMap;

    #[test]
    fn test_casher_basically_works() {
        let mut word_len = Casher::new(|word: &str| word.len());
        let hello = word_len.on("hello");

        // Test function returns correctly
        assert_eq!(hello, 5);

        // Test HashMap is correct length
        assert_eq!(word_len.values.len(), 1);

        // Test HashMap has correct value after one insert
        let mut test_map = HashMap::new();
        test_map.insert("hello", 5);
        assert_eq!(word_len.values, test_map);

        // Test HashMap has correct value after duplicate insert
        word_len.on("hello");
        assert_eq!(word_len.values, test_map);

        // Test HashMap has correct values after unique input
        word_len.on("wazzup");
        test_map.insert("wazzup", 6);
        assert_eq!(word_len.values, test_map);
    }

    #[test]
    fn testing_call_and_replace_and_forget() {
        // I'm honestly not sure how to write a good test
        // for call_and_replace, but if you have an idea let
        // me know and I'll put it here and in the docs

        // Testing forget
        let mut func = Casher::new(|x: usize| x * 2);
        func.on(27);
        func.on(54);
        func.on(222);
        func.forget();

        assert_eq!(func.values, HashMap::new());
    }
}
