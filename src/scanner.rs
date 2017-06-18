// Generic utility for reading data from standard input, based on
// http://codeforces.com/contest/702/submission/19589375
use std::io;

pub struct Scanner<B> {
    buffer: Vec<String>,
    reader: B,
}

impl<B: io::BufRead> Scanner<B> {
    pub fn new(reader: B) -> Self {
        Self {
            buffer: Vec::new(),
            reader: reader,
        }
    }

    // Use "turbofish" syntax next::<T>() to select data type of next token.
    pub fn next<T: ::std::str::FromStr>(&mut self) -> T
    where
        T::Err: ::std::fmt::Debug,
    {
        if let Some(front) = self.buffer.pop() {
            front.parse::<T>().expect(&front)
        } else {
            let mut input = String::new();
            self.reader.read_line(&mut input).expect("Line not read");
            self.buffer = input.split_whitespace().rev().map(String::from).collect();
            self.next()
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
        let x = scan.next::<i32>();
        let y = scan.next::<i32>();
        assert_eq!(x - y, 42);
    }

    #[test]
    fn test_stdin() {
        let stdin = io::stdin();
        let mut scan = Scanner::new(stdin.lock());
        if false {
            let _ = scan.next::<i32>();
        }
    }

    #[test]
    #[should_panic]
    fn test_file() {
        let file = ::std::fs::File::open("asdf.txt").expect("File not found");
        let mut scan = Scanner::new(io::BufReader::new(file));
        let _ = scan.next::<i32>();
    }
}
