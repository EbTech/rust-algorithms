// Palindrome substrings in O(n), Manacher's algorithm
// length of odd palin centred at s[i] is len[2*i]
// even palin btwn s[i],s[i+1]: len[2*i+1]
// TODO: check for underflows
// Alternative version:
// for c in (..)
//   while bla
//     len[c] += 2;
//     if len[c]-r+i == len[c-r+i] { len.push(len[c]-r+i); r += 1; }
fn find_pals(text: &[u8]) -> Vec<usize> {
    let mut len = Vec::with_capacity(2*text.len() - 1); 
    len.push(1); len.push(0);
    let mut i = 1;
    while i < 2*text.len() - 2 {
        let max_len = ::std::cmp::min(i+1, 2*text.len()-1-i);
        while len[i] < max_len && text[(i-len[i]-1)/2] == text[(i+len[i]+1)/2] {
            len[i] += 2;
        }
        let mut d = 1;
        while len[i-d] < len[i]-d { len[i+d] = len[i-d]; d += 1; }
        len[i+d] = len[i]-d;
        i += d;
    }
    len
}

// fail[i] = len of longest proper prefix-suffix of pat[0...i]
fn kmp_init(pat: &[u8]) -> Vec<usize> {
    let mut fail = Vec::with_capacity(pat.len());
    fail.push(0);
    let mut len = 0;
    for &ch in &pat[1..] {
        while len > 0 && pat[len] != ch { len = fail[len-1]; }
        if pat[len] == ch { len += 1; }
        fail.push(len);
    }
    fail
}

// matches[i] = len of longest prefix of pat matching with suffix of text[0...i]
fn kmp_match(text: &[u8], pat: &[u8]) -> Vec<usize> {
    let fail = kmp_init(pat);
    let mut matches = Vec::with_capacity(text.len());
    let mut len = 0;
    for &ch in text {
        if len == pat.len() { len = fail[len-1]; }
        while len > 0 && pat[len] != ch { len = fail[len-1]; }
        if pat[len] == ch { len += 1; }
        matches.push(len);
    }
    matches
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test_string() {
        let text = "abcbc".as_bytes();
        let pat = "bc".as_bytes();
        let matches = kmp_match(text, pat);
        //let pal_len = find_pals(text);
        assert_eq!(matches, vec![0, 1, 2, 1, 2]);
        //assert_eq!(pal_len, vec![1, 0, 1, 0, 3, 0, 3, 0, 1]);
    }
}
