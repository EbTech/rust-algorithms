//! String processing algorithms.

/// Data structure for Knuth-Morris-Pratt string matching against a pattern.
pub struct Matcher<'a> {
    /// The string pattern to search for.
    pub pattern: &'a [u8],
    /// KMP match failure automaton. fail[i] is the length of the longest
    /// proper prefix-suffix of pattern[0...i].
    pub fail: Vec<usize>,
}

impl<'a> Matcher<'a> {
    /// Precomputes the automaton that allows linear-time string matching.
    ///
    /// # Panics
    ///
    /// Panics if pattern is empty.
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
        Self { pattern, fail }
    }

    /// KMP algorithm, sets matches[i] = length of longest prefix of pattern
    /// matching a suffix of text[0...i].
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

/// Suffix array data structure, useful for a variety of string queries.
pub struct SuffixArray {
    /// The suffix array itself, holding suffix indices in sorted order.
    pub sfx: Vec<usize>,
    /// rank[i][j] = rank of the j'th suffix, considering only 2^i chars.
    pub rank: Vec<Vec<usize>>,
}

impl SuffixArray {
    /// O(n + max_key) stable sort on an input that is a permutation of (0..n).
    fn counting_sort<I>(p_gen: I, keys: &[usize], max_key: usize) -> Vec<usize>
    where
        I: DoubleEndedIterator<Item = usize>,
    {
        let mut counts = vec![0; max_key];
        for &k in keys {
            counts[k] += 1;
        }
        let mut total = 0;
        for c in counts.iter_mut() {
            total += *c;
            *c = total;
        }
        let mut result = vec![0; total];
        for p in p_gen.rev() {
            let c = &mut counts[keys[p]];
            *c -= 1;
            result[*c] = p;
        }
        result
    }

    /// Suffix array construction in O(n log n) time. Makes some unnecessary Vec clones
    /// and initializations, so there's room to optimize.
    pub fn new(text: &[u8]) -> Self {
        let n = text.len();
        let mut rank = vec![text.into_iter().map(|&ch| ch as usize).collect::<Vec<_>>()];
        let mut sfx = Self::counting_sort(0..n, rank.last().unwrap(), 256);
        // Invariant at the start of every loop iteration:
        // suffixes are sorted according to the first skip characters.
        for skip in (0..).map(|i| 1 << i).take_while(|&skip| skip < n) {
            let prev_rank = rank.last().unwrap().clone();
            let mut cur_rank = prev_rank.clone();

            let p_gen = (n - skip..n).chain(sfx.into_iter().filter_map(|p| p.checked_sub(skip)));
            sfx = Self::counting_sort(p_gen, &prev_rank, n.max(256));

            let mut prev = sfx[0];
            cur_rank[prev] = 0;
            for &p in sfx.iter().skip(1) {
                if prev.max(p) + skip < n && prev_rank[prev] == prev_rank[p] &&
                    prev_rank[prev + skip] == prev_rank[p + skip]
                {
                    cur_rank[p] = cur_rank[prev];
                } else {
                    cur_rank[p] = cur_rank[prev] + 1;
                }
                prev = p;
            }
            rank.push(cur_rank);
        }
        Self { sfx, rank }
    }

    /// Computes the length of longest common prefix of text[i..] and text[j..].
    pub fn longest_common_prefix(&self, mut i: usize, mut j: usize) -> usize {
        let mut len = 0;
        for (k, rank) in self.rank.iter().enumerate().rev() {
            if rank[i] == rank[j] {
                i += 1 << k;
                j += 1 << k;
                len += 1 << k;
                if i.max(j) >= self.sfx.len() {
                    break;
                }
            }
        }
        len
    }
}

/// Manacher's algorithm for computing palindrome substrings in linear time.
/// pal[2*i] = odd length of palindrome centred at text[i].
/// pal[2*i+1] = even length of palindrome centred at text[i+0.5].
///
/// # Panics
///
/// Panics if text is empty.
pub fn palindromes(text: &[u8]) -> Vec<usize> {
    let mut pal = Vec::with_capacity(2 * text.len() - 1); // only mutable var!
    pal.push(1);
    while pal.len() < pal.capacity() {
        let i = pal.len() - 1;
        let max_len = (i + 1).min(pal.capacity() - i);
        while pal[i] < max_len && text[(i - pal[i] - 1) / 2] == text[(i + pal[i] + 1) / 2] {
            pal[i] += 2;
        }
        if pal[i] < 2 {
            let a = 1 - pal[i];
            pal.push(a);
        } else {
            for d in 1.. {
                let (a, b) = (pal[i - d], pal[i] - d);
                if a < b {
                    pal.push(a);
                } else {
                    pal.push(b);
                    break;
                }
            }
        }
    }
    pal
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_kmp() {
        let text = "banana".as_bytes();
        let pattern = "ana".as_bytes();

        let matches = Matcher::new(pattern).kmp_match(text);

        assert_eq!(matches, vec![0, 1, 2, 3, 2, 3]);
    }

    #[test]
    fn test_suffix_array() {
        let text1 = "bobocel".as_bytes();
        let text2 = "banana".as_bytes();

        let sfx1 = SuffixArray::new(text1);
        let sfx2 = SuffixArray::new(text2);

        assert_eq!(sfx1.sfx, vec![0, 2, 4, 5, 6, 1, 3]);
        assert_eq!(sfx2.sfx, vec![5, 3, 1, 0, 4, 2]);

        assert_eq!(sfx1.longest_common_prefix(0, 2), 2);
        assert_eq!(sfx2.longest_common_prefix(1, 3), 3);

        // Check that sfx and rank.last() are essentially inverses of each other.
        for (p, &r) in sfx1.rank.last().unwrap().iter().enumerate() {
            assert_eq!(sfx1.sfx[r], p);
        }
        for (p, &r) in sfx2.rank.last().unwrap().iter().enumerate() {
            assert_eq!(sfx2.sfx[r], p);
        }
    }

    #[test]
    fn test_palindrome() {
        let text = "banana".as_bytes();

        let pal_len = palindromes(text);

        assert_eq!(pal_len, vec![1, 0, 1, 0, 3, 0, 5, 0, 3, 0, 1]);
    }
}
