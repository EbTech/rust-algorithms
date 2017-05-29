

pub struct Scanner {
    buffer: ::std::collections::VecDeque<String>
}

impl Scanner {
    pub fn new() -> Scanner {
        Scanner {
            buffer: ::std::collections::VecDeque::new()
        }
    }

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

