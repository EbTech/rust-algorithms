# Algorithm Cookbook in Rust

A collection of classic data structures and algorithms, emphasizing beauty and clarity over full generality. As such, this should not be viewed as a blackbox *library*, but as a whitebox *cookbook* demonstrating the translation of abstract concepts into executable code. I hope it will be useful to students and educators, as well as competition programmers.

This repository is distributed under the [MIT License](LICENSE). The license text need not be included in contest submissions, though I would appreciate linking back to this repo for others to find. Enjoy!

## For Students and Educators

When learning a new algorithm or data structure, it's often helpful to see or play with a concrete implementation. As such, this repository catalogues several classic algorithms in their simplest forms.

In addition, the Rust language has outstanding pedagogical attributes. Its compiler acts as a teacher, enforcing strict discipline while pointing to clearer ways to structure one's logic.

## For Competition Programmers

The original intent of this project was to build a reference for use in programming competitions such as [Codeforces](http://codeforces.com) and the [Google Code Jam](https://code.google.com/codejam). As a result, it contains algorithms that are frequently useful to have in one's toolkit, with an emphasis on making the code concise and easy to modify under time pressure.

Most competitive programmers use C/C++ because it allows for fast coding as well as fast execution. However, these languages are notoriously unsafe, wasting a considerable share of the contestant's time and attention on accident prevention and debugging. Java is the next most popular choice, offering a bit of safety at some expense to coding and execution speed. To my delight, I found that Rust provides a lot more safety than Java without the visual clutter, and it's *fast*. A proficient Rust programmer stands to gain a competitive advantage as well as a more pleasant experience!

## Programming Language Advocacy

My other goal is to appeal to developers who feel, as I once did, trapped between the lesser of headaches (e.g., C++ and Java), to raise awareness that *it doesn't have to be this way*. Rather than trying to persuade you with words, this repository aims to show by example and ease the learning curve a bit. See [Jim Blandy's *Why Rust?*](http://www.oreilly.com/programming/free/files/why-rust.pdf) for a brief introduction, or just [dive in](https://www.rust-lang.org)!

## Contents

- [Basic graph representations](src/graph/mod.rs): adjacency lists, minimum spanning tree, Euler path, disjoint set union 
- [Network flows](src/graph/flow.rs): Dinic's blocking flow, Hopcroft-Karp bipartite matching, min cost max flow
- [Connected components](src/graph/connectivity.rs): 2-edge-, 2-vertex- and strongly connected components, bridges, articulation points, topological sort, 2-SAT
- [Associative range query](src/arq_tree.rs): known colloquially as *segtrees*
- [Math](src/math.rs): Euclid's GCD algorithm, Bezout's identity
- [Scanner](src/scanner.rs): utility for reading input data
- [String processing](src/string_proc.rs): Knuth-Morris-Pratt string matching, Manacher's palindrome search
