# Contest Algorithms in Rust

[![Crates.io Version](https://img.shields.io/crates/v/contest-algorithms.svg)](https://crates.io/crates/contest-algorithms)
[![Documentation](https://docs.rs/contest-algorithms/badge.svg)](https://docs.rs/contest-algorithms)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/bevyengine/bevy/blob/master/LICENSE)
[![Crates.io Downloads](https://img.shields.io/crates/d/contest-algorithms.svg)](https://crates.io/crates/contest-algorithms)
[![Build Status](https://travis-ci.org/EbTech/rust-algorithms.svg?branch=master)](https://travis-ci.org/EbTech/rust-algorithms)
[![Gitter](https://badges.gitter.im/rust-algos/community.svg)](https://gitter.im/rust-algos/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

A collection of classic data structures and algorithms, emphasizing usability, beauty and clarity over full generality. As such, this should be viewed not as a blackbox *library*, but as a whitebox *cookbook* demonstrating the design and implementation of algorithms. I hope it will be useful to students and educators, as well as fans of algorithmic programming contests.

This repository is distributed under the [MIT License](LICENSE). Contest submissions need not include the license text. Enjoy!

## For Students and Educators

When learning a new algorithm or data structure, it's often helpful to see or play with a concrete implementation. As such, this repository catalogues several classic algorithms in their simplest forms.

In addition, the Rust language has outstanding pedagogical attributes. Its compiler acts as a teacher, enforcing strict discipline while pointing to clearer ways to structure one's logic.

## For Programming Contests

The original intent of this project was to build a reference for use in programming contests. As a result, it contains algorithms that are frequently useful to have in one's toolkit, with an emphasis on code that is concise and easy to modify under time pressure.

Most competitive programmers rely on C++ for its fast execution time. However, it's notoriously unsafe, diverting a considerable share of the contestant's time and attention on mistake prevention and debugging. Java is the next most popular choice, offering a little safety at some expense to speed of coding and execution.

To my delight, I found that Rust eliminates entire classes of bugs, while reducing visual clutter to make the rest easier to spot. And, it's *fast*. There's a learning curve, to be sure. However, a proficient Rust programmer stands to gain a competitive advantage as well as a more pleasant experience!

Some contest sites and online judges that support Rust:
- [Codeforces](https://codeforces.com)
- [AtCoder](https://atcoder.jp)
- [Kattis](https://open.kattis.com/help/rust)
- [SPOJ](https://www.spoj.com/)
- [LeetCode](https://leetcode.com/contest)
- [HackerRank](https://www.hackerrank.com/contests)
- [Timus](http://acm.timus.ru/help.aspx?topic=rust)

The following support pre-2018 versions of Rust:
- [Google Kick Start and Code Jam](https://codingcompetitions.withgoogle.com)

For help in getting started, you may check out [some of my past submissions](https://codeforces.com/contest/1168/submission/55200038).

## Programming Language Advocacy

My other goal is to appeal to developers who feel limited by ancient (yet still mainstream) programming languages, by demonstrating the power of modern techniques.

Rather than try to persuade you with words, this repository aims to show by example. If you'd like to learn the language, I recommend [the official book](https://doc.rust-lang.org/book/) or [Programming Rust](https://www.amazon.com/Programming-Rust-Fast-Systems-Development-dp-1492052590/dp/1492052590).

# Contents

## [Graphs](src/graph/)

### [Graph representations](src/graph/mod.rs)

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
- Hopcroft-Karp bipartite matching
- Minimum cost maximum flow

## [Math](src/math/)

### [Number theory](src/math/mod.rs)

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

