# Algorithms Cookbook in Rust

A collection of famous data structures and algorithms, emphasizing beauty and clarity over full generality. As such, this should not be viewed as a blackbox *library*, but as a whitebox *cookbook* demonstrating how these abstract concepts are translated into [Rust](https://www.rust-lang.org) code. I hope it will be useful to students and educators, as well as competition programmers.

This repository is distributed under the [MIT License](LICENSE). The license text need not be included in contest submissions, though I would appreciate linking back to this repo for others to find. Enjoy!

## Academic Programming Competitions
The original intent of this project was to build a reference for use in programming competitions such as [Codeforces](http://codeforces.com) and the [Google Code Jam](https://code.google.com/codejam). As a result, it contains algorithms that are frequently useful to have in one's toolkit, with an emphasis on making the code concise and easy to modify in time-constrained settings.

Most competitive programmers use C/C++ because it allows for fast coding as well as fast execution. However, these languages are notoriously unsafe, wasting a considerable share of the contestant's time and attention on accident prevention. Java is the next most popular choice, offering a bit of safety at some expense to coding and execution speed. To my delight, I found that Rust provides a lot more safety than Java without the visual clutter, and it's *fast*. 

## Programming Language Advocacy

My other goal is to show the world that C++ kinda sucks, and that *it doesn't have to be that way*. Rather than trying to persuade you with words, this repository aims to show by example.

## Contents

- [Basic graph representations](src/graph/mod.rs): adjacency lists, minimum spanning tree, Euler path, disjoint set union 
- [Network flows](src/graph/flow.rs): Dinic's blocking flow, min cost max flow, bipartite matching
- [Connected components](src/graph/connectivity.rs): 2-edge-, 2-vertex- and strongly connected components, topological sort, 2-SAT
- [Associative range query](src/arq_tree.rs): known colloquially as *segtrees*
- [Math](src/math.rs): Euclid's GCD algorithm, Bezout's identity
- [Scanner](src/scanner.rs): utility for reading data from standard input
- [String processing](src/string_proc.rs): Knuth-Morris-Pratt string matching, Manacher's palindrome search
