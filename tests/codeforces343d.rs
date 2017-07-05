//! Solves [Water Tree](http://codeforces.com/contest/343/problem/D).
//! To make a self-contained file for contest submission, dump each desired
//! module's contents directly here instead of the use statements.
//! Also, replace io::Cursor with io::stdin as shown in scanner.rs.
extern crate algorithms;
use algorithms::arq_tree::{ArqTree, AssignAdd};
use algorithms::graph::Graph;
use algorithms::scanner::Scanner;

const SAMPLE_INPUT: &str = "\
5
1 2
5 1
2 3
4 2
12
1 1
2 3
3 1
3 2
3 3
3 4
1 2
2 4
3 1
3 3
3 4
3 5
";
const SAMPLE_OUTPUT: &str = "\
0
0
0
1
0
1
0
1
";

fn dfs(
    graph: &Graph,
    u: usize,
    l: &mut [usize],
    r: &mut [usize],
    p: &mut [usize],
    time: &mut usize,
) {
    *time += 1;
    l[u] = *time;

    for (_, v) in graph.adj_list(u) {
        if l[v] == 0 {
            p[v] = l[u];
            dfs(graph, v, l, r, p, time);
        }
    }

    r[u] = *time;
}

fn main1() {
    let cursor = std::io::Cursor::new(SAMPLE_INPUT);
    let mut scan = Scanner::new(cursor);
    let mut out = String::new();

    let n = scan.next::<usize>();
    let mut tree = Graph::new(n, 2 * (n - 1));
    for _ in 1..n {
        let u = scan.next::<usize>() - 1;
        let v = scan.next::<usize>() - 1;
        tree.add_undirected_edge(u, v);
    }

    let mut l = vec![0; n];
    let mut r = vec![0; n];
    let mut p = vec![0; n];
    dfs(&tree, 0, &mut l, &mut r, &mut p, &mut 0);

    let mut arq = ArqTree::<AssignAdd>::new(vec![(0, 1); n + 1]);
    let q = scan.next::<usize>();
    for _ in 0..q {
        let c = scan.next::<usize>();
        let v = scan.next::<usize>() - 1;
        let (sum, len) = arq.query(l[v], r[v]);
        if c == 1 {
            if sum != len {
                arq.modify(p[v], p[v], &0);
                arq.modify(l[v], r[v], &1);
            }
        } else if c == 2 {
            arq.modify(l[v], l[v], &0);
        } else {
            use std::fmt::Write;
            writeln!(&mut out, "{}", if sum == len { 1 } else { 0 }).unwrap();
        }
    }

    assert_eq!(out, SAMPLE_OUTPUT);
}

#[test]
fn main() {
    // If your contest solution requires a lot of stack space, make sure to
    // run it in a custom thread like this.
    std::thread::Builder::new()
        .stack_size(50_000_000)
        .spawn(main1)
        .unwrap()
        .join()
        .unwrap();
}
