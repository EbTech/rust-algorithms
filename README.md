# Contest Algorithms in Rust

[![Build Status](https://travis-ci.org/EbTech/rust-algorithms.svg?branch=master)](https://travis-ci.org/EbTech/rust-algorithms)
[![Latest Version](https://img.shields.io/crates/v/contest-algorithms.svg)](https://crates.io/crates/contest-algorithms)

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
- [Kattis](https://open.kattis.com/help/rust)
- [SPOJ](https://www.spoj.com/)
- [LeetCode](https://leetcode.com/contest)
- [HackerRank](https://www.hackerrank.com/contests)

The following support pre-2018 versions of Rust:
- [Google Kick Start and Code Jam](https://codingcompetitions.withgoogle.com)
- [AtCoder](https://atcoder.jp)
- [Timus](http://acm.timus.ru/help.aspx?topic=rust)

For help in getting started, you may check out [some of my past submissions](https://codeforces.com/contest/1168/submission/55200038).

## Programming Language Advocacy

My other goal is to appeal to developers who feel limited by ancient (yet still mainstream) programming languages, by demonstrating the power of modern techniques.

Rather than try to persuade you with words, this repository aims to show by example. If you're new to Rust, see [Jim Blandy's *Why Rust?*](http://www.oreilly.com/programming/free/files/why-rust.pdf) for a brief introduction, or just [dive in!](https://doc.rust-lang.org/book/)

## Contents

- [Basic graph representations](src/graph/mod.rs): adjacency lists, disjoint set union
- [Elementary graph algorithms](src/graph/util.rs): minimum spanning tree, Euler path, Dijkstra's algorithm, DFS iteration
- [Connected components](src/graph/connectivity.rs): 2-edge-, 2-vertex- and strongly connected components, bridges, articulation points, topological sort, 2-SAT
- [Network flows](src/graph/flow.rs): Dinic's blocking flow, Hopcroft-Karp bipartite matching, min cost max flow
- [Number theory](src/math/mod.rs): canonical solution to Bezout's identity, Miller's primality test
- [FFT](src/math/fft.rs): fast Fourier transform, number theoretic transform, convolution
- [Arithmetic](src/math/num.rs): rational and complex numbers, linear algebra, safe modular arithmetic
- [Ordering algorithms](src/order.rs): binary search, mergesort, coordinate compression, online convex hull trick
- [Associative range query](src/range_query): static and dynamic ARQ trees (a.k.a. segtrees), Mo's query square root decomposition
- [Scanner](src/scanner.rs): utility for reading input data ergonomically
- [String processing](src/string_proc.rs): Knuth-Morris-Pratt and Aho-Corasick string matching, suffix array, Manacher's linear-time palindrome search

