// Generic utility for reading data from standard input.
pub struct Scanner {
    buffer: ::std::collections::VecDeque<String>
}

impl Scanner {
    pub fn new() -> Scanner {
        Scanner {
            buffer: ::std::collections::VecDeque::new()
        }
    }
    
    // Use "turbofish" syntax next::<T>() to select data type of next token.
    pub fn next<T: ::std::str::FromStr>(&mut self) -> T {
        while self.buffer.is_empty() {
            let mut input = String::new();
            ::std::io::stdin().read_line(&mut input).ok();
            self.buffer = input.split_whitespace().map(String::from).collect();
        }
        let front = self.buffer.pop_front().unwrap();
        front.parse::<T>().ok().unwrap()
    }
}

