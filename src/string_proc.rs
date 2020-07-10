//! String processing algorithms.
use std::cmp::{max, min};
use std::collections::{hash_map::Entry, HashMap, VecDeque};

/// Prefix trie, easily augmentable by adding more fields and/or methods
pub struct Trie<C: std::hash::Hash + Eq> {
    links: Vec<HashMap<C, usize>>,
}

impl<C: std::hash::Hash + Eq> Default for Trie<C> {
    /// Creates an empty trie with a root node.
    fn default() -> Self {
        Self {
            links: vec![HashMap::new()],
        }
    }
}

impl<C: std::hash::Hash + Eq> Trie<C> {
    /// Inserts a word into the trie, and returns the index of its node.
    pub fn insert(&mut self, word: impl IntoIterator<Item = C>) -> usize {
        let mut node = 0;

        for ch in word {
            let len = self.links.len();
            node = match self.links[node].entry(ch) {
                Entry::Occupied(entry) => *entry.get(),
                Entry::Vacant(entry) => {
                    entry.insert(len);
                    self.links.push(HashMap::new());
                    len
                }
            }
        }
        node
    }

    /// Finds a word in the trie, and returns the index of its node.
    pub fn get(&self, word: impl IntoIterator<Item = C>) -> Option<usize> {
        let mut node = 0;
        for ch in word {
            node = *self.links[node].get(&ch)?;
        }
        Some(node)
    }
}

/// Single-pattern matching with the Knuth-Morris-Pratt algorithm
pub struct Matcher<'a, C: Eq> {
    /// The string pattern to search for.
    pub pattern: &'a [C],
    /// KMP match failure automaton. fail[i] is the length of the longest
    /// proper prefix-suffix of pattern[0..=i].
    pub fail: Vec<usize>,
}

impl<'a, C: Eq> Matcher<'a, C> {
    /// Precomputes the automaton that allows linear-time string matching.
    ///
    /// # Example
    ///
    /// ```
    /// use contest_algorithms::string_proc::Matcher;
    /// let utf8_string = "hello";
    ///
    /// let match_from_byte_literal = Matcher::new(b"hello");
    ///
    /// let match_from_bytes = Matcher::new(utf8_string.as_bytes());
    ///
    /// let vec_char: Vec<char> = utf8_string.chars().collect();
    /// let match_from_chars = Matcher::new(&vec_char);
    ///
    /// let vec_int = vec![4, -3, 1];
    /// let match_from_ints = Matcher::new(&vec_int);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if pattern is empty.
    pub fn new(pattern: &'a [C]) -> Self {
        let mut fail = Vec::with_capacity(pattern.len());
        fail.push(0);
        let mut len = 0;
        for ch in &pattern[1..] {
            while len > 0 && pattern[len] != *ch {
                len = fail[len - 1];
            }
            if pattern[len] == *ch {
                len += 1;
            }
            fail.push(len);
        }
        Self { pattern, fail }
    }

    /// KMP algorithm, sets match_lens[i] = length of longest prefix of pattern
    /// matching a suffix of text[0..=i].
    pub fn kmp_match(&self, text: &[C]) -> Vec<usize> {
        let mut match_lens = Vec::with_capacity(text.len());
        let mut len = 0;
        for ch in text {
            if len == self.pattern.len() {
                len = self.fail[len - 1];
            }
            while len > 0 && self.pattern[len] != *ch {
                len = self.fail[len - 1];
            }
            if self.pattern[len] == *ch {
                len += 1;
            }
            match_lens.push(len);
        }
        match_lens
    }
}

/// Multi-pattern matching with the Aho-Corasick algorithm
pub struct MultiMatcher<C: std::hash::Hash + Eq> {
    /// A prefix trie storing the string patterns to search for.
    pub trie: Trie<C>,
    /// Stores which completed pattern string each node corresponds to.
    pub pat_id: Vec<Option<usize>>,
    /// Aho-Corasick failure automaton. fail[i] is the node corresponding to the
    /// longest prefix-suffix of the node corresponding to i.
    pub fail: Vec<usize>,
    /// Shortcut to the next match along the failure chain, or to the root.
    pub fast: Vec<usize>,
}

impl<C: std::hash::Hash + Eq> MultiMatcher<C> {
    fn next(trie: &Trie<C>, fail: &[usize], mut node: usize, ch: &C) -> usize {
        loop {
            if let Some(&child) = trie.links[node].get(ch) {
                return child;
            } else if node == 0 {
                return 0;
            }
            node = fail[node];
        }
    }

    /// Precomputes the automaton that allows linear-time string matching.
    /// If there are duplicate patterns, all but one copy will be ignored.
    pub fn new(patterns: Vec<impl IntoIterator<Item = C>>) -> Self {
        let mut trie = Trie::default();
        let pat_nodes: Vec<usize> = patterns.into_iter().map(|pat| trie.insert(pat)).collect();

        let mut pat_id = vec![None; trie.links.len()];
        for (i, node) in pat_nodes.into_iter().enumerate() {
            pat_id[node] = Some(i);
        }

        let mut fail = vec![0; trie.links.len()];
        let mut fast = vec![0; trie.links.len()];
        let mut q: VecDeque<usize> = trie.links[0].values().cloned().collect();

        while let Some(node) = q.pop_front() {
            for (ch, &child) in &trie.links[node] {
                let nx = Self::next(&trie, &fail, fail[node], &ch);
                fail[child] = nx;
                fast[child] = if pat_id[nx].is_some() { nx } else { fast[nx] };
                q.push_back(child);
            }
        }

        Self {
            trie,
            pat_id,
            fail,
            fast,
        }
    }

    /// Aho-Corasick algorithm, sets match_nodes[i] = node corresponding to
    /// longest prefix of some pattern matching a suffix of text[0..=i].
    pub fn ac_match(&self, text: &[C]) -> Vec<usize> {
        let mut match_nodes = Vec::with_capacity(text.len());
        let mut node = 0;
        for ch in text {
            node = Self::next(&self.trie, &self.fail, node, &ch);
            match_nodes.push(node);
        }
        match_nodes
    }

    /// For each non-empty match, returns where in the text it ends, and the index
    /// of the corresponding pattern.
    pub fn get_end_pos_and_pat_id(&self, match_nodes: &[usize]) -> Vec<(usize, usize)> {
        let mut res = vec![];
        for (text_pos, &(mut node)) in match_nodes.iter().enumerate() {
            while node != 0 {
                if let Some(id) = self.pat_id[node] {
                    res.push((text_pos + 1, id));
                }
                node = self.fast[node];
            }
        }
        res
    }
}

/// Suffix array data structure, useful for a variety of string queries.
pub struct SuffixArray {
    /// The suffix array itself, holding suffix indices in sorted order.
    pub sfx: Vec<usize>,
    /// rank[i][j] = rank of the j'th suffix, considering only 2^i chars.
    /// In other words, rank[i] is a ranking of the substrings text[j..j+2^i].
    pub rank: Vec<Vec<usize>>,
}

impl SuffixArray {
    /// O(n + max_key) stable sort on the items generated by vals.
    /// Items v in vals are sorted according to val_to_key[v].
    fn counting_sort(
        vals: impl Iterator<Item = usize> + Clone,
        val_to_key: &[usize],
        max_key: usize,
    ) -> Vec<usize> {
        let mut counts = vec![0; max_key];
        for v in vals.clone() {
            counts[val_to_key[v]] += 1;
        }
        let mut total = 0;
        for c in counts.iter_mut() {
            total += *c;
            *c = total - *c;
        }
        let mut result = vec![0; total];
        for v in vals {
            let c = &mut counts[val_to_key[v]];
            result[*c] = v;
            *c += 1;
        }
        result
    }

    /// Suffix array construction in O(n log n) time.
    pub fn new(text: &[u8]) -> Self {
        let n = text.len();
        let init_rank = text.iter().map(|&ch| ch as usize).collect::<Vec<_>>();
        let mut sfx = Self::counting_sort(0..n, &init_rank, 256);
        let mut rank = vec![init_rank];
        // Invariant at the start of every loop iteration:
        // suffixes are sorted according to the first skip characters.
        for skip in (0..).map(|i| 1 << i).take_while(|&skip| skip < n) {
            let prev_rank = rank.last().unwrap();
            let mut cur_rank = prev_rank.clone();

            let pos = (n - skip..n).chain(sfx.into_iter().filter_map(|p| p.checked_sub(skip)));
            sfx = Self::counting_sort(pos, &prev_rank, max(n, 256));

            let mut prev = sfx[0];
            cur_rank[prev] = 0;
            for &cur in sfx.iter().skip(1) {
                if max(prev, cur) + skip < n
                    && prev_rank[prev] == prev_rank[cur]
                    && prev_rank[prev + skip] == prev_rank[cur + skip]
                {
                    cur_rank[cur] = cur_rank[prev];
                } else {
                    cur_rank[cur] = cur_rank[prev] + 1;
                }
                prev = cur;
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
                if max(i, j) >= self.sfx.len() {
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
pub fn palindromes<T: Eq>(text: &[T]) -> Vec<usize> {
    let mut pal = Vec::with_capacity(2 * text.len() - 1);
    pal.push(1);
    while pal.len() < pal.capacity() {
        let i = pal.len() - 1;
        let max_len = min(i + 1, pal.capacity() - i);
        while pal[i] < max_len && text[(i - pal[i] - 1) / 2] == text[(i + pal[i] + 1) / 2] {
            pal[i] += 2;
        }
        if let Some(a) = 1usize.checked_sub(pal[i]) {
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
    fn test_kmp_matching() {
        let text = b"banana";
        let pattern = b"ana";

        let matches = Matcher::new(pattern).kmp_match(text);

        assert_eq!(matches, vec![0, 1, 2, 3, 2, 3]);
    }

    #[test]
    fn test_ac_matching() {
        let text = b"banana bans, apple benefits.";
        let dict = vec![
            "banana".bytes(),
            "benefit".bytes(),
            "banapple".bytes(),
            "ban".bytes(),
            "fit".bytes(),
        ];

        let matcher = MultiMatcher::new(dict);
        let match_nodes = matcher.ac_match(text);
        let end_pos_and_id = matcher.get_end_pos_and_pat_id(&match_nodes);

        assert_eq!(
            end_pos_and_id,
            vec![(3, 3), (6, 0), (10, 3), (26, 1), (26, 4)]
        );
    }

    #[test]
    fn test_suffix_array() {
        let text1 = b"bobocel";
        let text2 = b"banana";

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
        let text = b"banana";

        let pal_len = palindromes(text);

        assert_eq!(pal_len, vec![1, 0, 1, 0, 3, 0, 5, 0, 3, 0, 1]);
    }
}
