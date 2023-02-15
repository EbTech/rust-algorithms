# Rust /\r5 /\lgorithmic/\

This is a collection of classic data structures and interesting algorithms, emphasizing clarity, elegance and understanding over generality and speed. One design criterion is to transparency: explain the names and concepts, annotate with complexity, be obvious in intention. Another is simplicity: the Rust ecosystem is full of well meaning crates, and I have intentionally decided to stick with data structures and algorithms only found in std for universality.

This is a fork of EbTech's amazing repo. My intention is to change the fundamental graph base classes to more in keeping with our intuition and understanding of graphs, do more data representation encapsulation, decouple algorithm implementation from knowledge of the internals of their datastructures, put in more tests of algorithms validity and speed,and when happy with that, start adding far more algorithms. 

My hope is that someday my kids will use this to learn about algorithms, something that I have always been obsessed about. Its also intended for students/teachers of algorithmica. I also think that its useful for competive programming; and that Rust is in many ways a good language for that. Except for the fact that its not easy to make a linked list...

Rust is a language that makes algorithms smaller and simpler, eliminating some deep bugs that are inevitable from the complexity in other languages (like c++). Its functional nature enforces elegance and compactness, its compiler assists in correctness.


![llama](https://user-images.githubusercontent.com/9121210/218507152-5a9646d5-c8bb-4937-acfb-8834410975fd.jpg)

Some contest sites and online judges that support Rust:
- [Codeforces](https://codeforces.com)
- [AtCoder](https://atcoder.jp)
- [Kattis](https://open.kattis.com/help/rust)
- [SPOJ](https://www.spoj.com/)
- [LeetCode](https://leetcode.com/contest)
- [HackerRank](https://www.hackerrank.com/contests)
- [Timus](http://acm.timus.ru/help.aspx?topic=rust)

# Contents

## [Graphs](src/graph/)

### [Graph representations](src/graph/graph.rs)

- Integer index-based adjacency list representation
- Disjoint set union

### [Elementary graph algorithms](src/graph/util.rs)

- Euler path and tour
- Kruskal's minimum spanning tree 
- Dijkstra's single-source shortest paths
- DFS pre-order traversal

### [Connected components](src/graph/connectivity.rs)

- Connected components
- Strongly connected components
- Bridges and 2-edge-connected components
- Articulation points and 2-vertex-connected components
- Topological sort
- 2-SAT solver

### [Network flows](src/graph/flow.rs)

- Dinic's blocking maximum flow
- Minimum cut
- Hopcroft-Karp bipartite maximum matching O(\sqrt(V)*E)
- Minimum cost maximum flow

## [Math](src/math/)

### [Number theory](src/math/division.rs)

- Greatest common divisor
- Canonical solution to Bezout's identity
- Miller's primality test

### [Generic FFT](src/math/fft.rs)

- Fast Fourier transform
- Number theoretic transform
- Convolution

### [Arithmetic](src/math/num.rs)

- Rational numbers
- Complex numbers
- Linear algebra
- Safe modular arithmetic

## [Ordering and search](src/order.rs)

- Comparator for `PartialOrd`
- Binary search: drop-in replacements for C++ `lower_bound()`/`upper_bound()`
- Merge and mergesort
- Coordinate compression
- Online convex hull trick (update and query the upper envelope of a set of lines)

## [Associative range query](src/range_query)

- Statically allocated binary indexed ARQ tree (a.k.a. generic segtree with lazy propagation)
- Dynamically allocated ARQ tree, optionally sparse and persistent
- Mo's algorithm (a.k.a. query square root decomposition)

## [Scanner](src/scanner.rs)

- Utility for reading input data ergonomically
- File and standard I/O examples

## [String processing](src/string_proc.rs)

- Generic trie
- Knuth-Morris-Pratt single-pattern string matching
- Aho-Corasick multi-pattern string matching
- Suffix array: O(n log n) construction using counting sort
- Longest common prefix
- Manacher's linear-time palindrome search

