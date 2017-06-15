// Data structure for Knuth-Morris-Pratt string matching against a pattern.
pub struct Matcher<'a> {
    pub pattern: &'a [u8],
    pub fail: Vec<usize>,
}

impl<'a> Matcher<'a> {
    // Sets fail[i] = length of longest proper prefix-suffix of pattern[0...i].
    pub fn new(pattern: &'a [u8]) -> Self {
        let mut fail = Vec::with_capacity(pattern.len());
        fail.push(0);
        let mut len = 0;
        for &ch in &pattern[1..] {
            while len > 0 && pattern[len] != ch {
                len = fail[len - 1];
            }
            if pattern[len] == ch {
                len += 1;
            }
            fail.push(len);
        }
        Self {
            pattern: pattern,
            fail: fail,
        }
    }

    // KMP algorithm, sets matches[i] = length of longest prefix of pattern
    // matching a suffix of text[0...i].
    pub fn kmp_match(&self, text: &[u8]) -> Vec<usize> {
        let mut matches = Vec::with_capacity(text.len());
        let mut len = 0;
        for &ch in text {
            if len == self.pattern.len() {
                len = self.fail[len - 1];
            }
            while len > 0 && self.pattern[len] != ch {
                len = self.fail[len - 1];
            }
            if self.pattern[len] == ch {
                len += 1;
            }
            matches.push(len);
        }
        matches
    }
}

// Manacher's algorithm for computing palindrome substrings in linear time.
// len[2*i] = odd length of palindrome centred at text[i].
// len[2*i+1] = even length of palindrome centred at text[i+0.5].
pub fn palindromes(text: &[u8]) -> Vec<usize> {
    let mut len = Vec::with_capacity(2 * text.len() - 1);
    len.push(1);
    while len.len() < len.capacity() {
        let i = len.len() - 1;
        let max_len = ::std::cmp::min(i + 1, len.capacity() - i);
        while len[i] < max_len && text[(i - len[i] - 1) / 2] == text[(i + len[i] + 1) / 2] {
            len[i] += 2;
        }
        if len[i] < 2 {
            let a = 1 - len[i];
            len.push(a);
        } else {
            for d in 1.. {
                let (a, b) = (len[i - d], len[i] - d);
                if a < b {
                    len.push(a);
                } else {
                    len.push(b);
                    break;
                }
            }
        }
    }
    len
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_string() {
        let text = "abcbc".as_bytes();
        let pattern = "bc".as_bytes();

        let matches = Matcher::new(pattern).kmp_match(text);
        assert_eq!(matches, vec![0, 1, 2, 1, 2]);

        let pal_len = palindromes(text);
        assert_eq!(pal_len, vec![1, 0, 1, 0, 3, 0, 3, 0, 1]);
    }
}
