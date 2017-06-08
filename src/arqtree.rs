// Associative Range Query Tree based on http://codeforces.com/blog/entry/18051
// Entries [0...size-1] are stored in t[size..2*size-1].
// The range operation must be associative: in this example, we use addition.
// In this example, the range operation assigns the value op to all entries.
pub struct ARQT {
    d: Vec<Option<i32>>,
    t: Vec<i32>,
    s: Vec<i32>
}

impl ARQT {
    pub fn new(size: usize) -> ARQT {
        let mut s = vec![1; 2*size];
        for i in (0..size).rev() {
            s[i] = s[i<<1] + s[i<<1|1];
        }
        ARQT {
            d: vec![None; size],
            t: vec![0; 2*size],                       // monoid identity
            s: s
        }
    }
    
    fn apply(&mut self, p: usize, op: i32) {
        self.t[p] = op * self.s[p];                   // hom application
        if p < self.d.len() { self.d[p] = Some(op); } // hom composition
    }
    
    fn push(&mut self, p: usize) {
        for s in (1..32).rev() {
            let i = p >> s;
            if let Some(op) = self.d[i] {
                self.apply(i<<1, op);
                self.apply(i<<1|1, op);
                self.d[i] = None;
            }
        }
    }
    
    fn pull(&mut self, mut p: usize) {
        while p > 1 {
            p >>= 1;
            if self.d[p] == None {
                self.t[p] = self.t[p<<1] + self.t[p<<1|1]; // monoid op
            }
        }
    }
    
    pub fn modify(&mut self, mut l: usize, mut r: usize, op: i32) {
        l += self.d.len(); r += self.d.len();
        let (l0, r0) = (l, r);
        self.push(l0); self.push(r0);
        while l <= r {
          if l & 1 == 1 { self.apply(l, op); l += 1; }
          if r & 1 == 0 { self.apply(r, op); r -= 1; }
          l >>= 1; r >>= 1;
        }
        self.pull(l0); self.pull(r0);
    }
    
    pub fn query(&mut self, mut l: usize, mut r: usize) -> i32 {
        l += self.d.len(); r += self.d.len();
        self.push(l); self.push(r);
        let mut res = 0;                                     // monoid identity
        while l <= r {
            if l & 1 == 1 { res = res + self.t[l]; l += 1; } // monoid op
            if r & 1 == 0 { res = self.t[r] + res; r -= 1; } // monoid op
            l >>= 1; r >>= 1;
        }
        res
    }
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test_arqt()
    {
        let mut arqt = ARQT::new(10);
        arqt.modify(1, 3, 10);
        arqt.modify(3, 5, 1);
        assert_eq!(arqt.query(0, 9), 23);
    }
}

