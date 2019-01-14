//! Generic utility for reading data from standard input, based on [voxl's
//! stdin wrapper](http://codeforces.com/contest/702/submission/19589375).
use std::io;

/// Reads white-space separated tokens one at a time.
pub struct Scanner<B> {
    reader: B,
    buffer: Vec<String>,
}

impl<B: io::BufRead> Scanner<B> {
    pub fn new(reader: B) -> Self {
        Self {
            reader,
            buffer: Vec::new(),
        }
    }

    /// Use "turbofish" syntax read::<T>() to select data type of next token.
    ///
    /// # Panics
    ///
    /// Panics if there's an I/O error or if the token cannot be parsed as T.
    pub fn read<T: ::std::str::FromStr>(&mut self) -> T
    where
        T::Err: ::std::fmt::Debug,
    {
        loop {
            if let Some(front) = self.buffer.pop() {
                return front.parse::<T>().expect(&front);
            }
            let mut input = String::new();
            self.reader.read_line(&mut input).expect("Line not read");
            self.buffer = input.split_whitespace().rev().map(String::from).collect();
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_fake_input() {
        let cursor = io::Cursor::new("44 2");
        let mut scan = Scanner::new(cursor);
        let x = scan.read::<i32>();
        let y = scan.read::<i32>();
        assert_eq!(x - y, 42);
    }

    #[test]
    fn test_stdin() {
        let stdin = io::stdin();
        let mut scan = Scanner::new(stdin.lock());
        if false {
            let _ = scan.read::<i32>();
        }
    }

    #[test]
    #[should_panic(expected = "File not found")]
    fn test_file() {
        let file = ::std::fs::File::open("asdf.txt").expect("File not found");
        let mut scan = Scanner::new(io::BufReader::new(file));
        let _ = scan.read::<i32>();
    }
}
