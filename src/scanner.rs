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

    /// Use "turbofish" syntax token::<T>() to select data type of next token.
    ///
    /// # Panics
    ///
    /// Panics if there's an I/O error or if the token cannot be parsed as T.
    pub fn token<T: std::str::FromStr>(&mut self) -> T {
        loop {
            if let Some(token) = self.buffer.pop() {
                return token.parse().ok().expect("Failed parse");
            }
            let mut input = String::new();
            self.reader.read_line(&mut input).expect("Failed read");
            self.buffer = input.split_whitespace().rev().map(String::from).collect();
        }
    }
}

pub fn scanner_from_file(filename: &str) -> Scanner<io::BufReader<std::fs::File>> {
    let file = std::fs::File::open(filename).expect("Input file not found");
    Scanner::new(io::BufReader::new(file))
}

pub fn writer_to_file(filename: &str) -> io::BufWriter<std::fs::File> {
    let file = std::fs::File::open(filename).expect("Output file not found");
    io::BufWriter::new(file)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_in_memory_io() {
        let cursor = io::Cursor::new("50 8");
        let mut scan = Scanner::new(cursor);
        let mut out = String::new();
        use std::fmt::Write; // needed for writeln!()

        let x = scan.token::<i32>();
        let y = scan.token::<i32>();
        writeln!(out, "Test {}", x - y).ok();

        assert_eq!(out, "Test 42\n");
    }

    #[test]
    fn test_compile_stdio() {
        let (stdin, stdout) = (io::stdin(), io::stdout());
        let mut scan = Scanner::new(stdin.lock());
        let mut out = io::BufWriter::new(stdout.lock());
        use io::Write; // needed for writeln!()

        if false {
            let x = scan.token::<i32>();
            let y = scan.token::<i32>();
            writeln!(out, "Test {}", x - y).ok();
        }
    }

    #[test]
    #[should_panic(expected = "Input file not found")]
    fn test_panic_file() {
        let mut scan = scanner_from_file("input_file.txt");
        let mut out = writer_to_file("output_file.txt");
        use io::Write; // needed for writeln!()

        let x = scan.token::<i32>();
        let y = scan.token::<i32>();
        writeln!(out, "Test {}", x - y).ok();
    }
}
