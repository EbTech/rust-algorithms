//! Solves [Water Tree](http://codeforces.com/contest/343/problem/D).
//! To make a self-contained file for contest submission, dump each desired
//! module's contents directly here instead of the use statements.
//! Also, replace io::Cursor with io::stdin as shown in scanner.rs.
extern crate contest_algorithms;
use contest_algorithms::arq_tree::{ArqTree, AssignSum};
use contest_algorithms::graph::Graph;
use contest_algorithms::scanner::Scanner;

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

#[test]
fn main() {
    use std::fmt::Write;
    let cursor = std::io::Cursor::new(SAMPLE_INPUT);
    let mut scan = Scanner::new(cursor);
    let out = &mut String::new();
    /* To read/write with stdin/stdout instead:
        use std::io::{BufWriter, stdin, stdout, Write};
        let stdin = stdin();
        let mut scan = Scanner::new(stdin.lock());
        let out = &mut BufWriter::new(stdout());
    */

    let n = scan.read::<usize>();
    let mut tree = Graph::new(n, 2 * (n - 1));
    for _ in 1..n {
        let u = scan.read::<usize>() - 1;
        let v = scan.read::<usize>() - 1;
        tree.add_undirected_edge(u, v);
    }

    let mut l = vec![0; n];
    let mut r = vec![0; n];
    let mut p = vec![0; n];
    dfs(&tree, 0, &mut l, &mut r, &mut p, &mut 0);

    let mut arq = ArqTree::<AssignSum>::new(vec![(0, 1); n + 1]);
    let q = scan.read::<usize>();
    for _ in 0..q {
        let c = scan.read::<usize>();
        let v = scan.read::<usize>() - 1;
        let (sum, len) = arq.query(l[v], r[v]);
        if c == 1 {
            if sum != len {
                arq.modify(p[v], p[v], &0);
                arq.modify(l[v], r[v], &1);
            }
        } else if c == 2 {
            arq.modify(l[v], l[v], &0);
        } else {
            let ans = if sum == len { 1 } else { 0 };
            writeln!(out, "{}", ans).ok();
        }
    }

    assert_eq!(out, SAMPLE_OUTPUT);
}
