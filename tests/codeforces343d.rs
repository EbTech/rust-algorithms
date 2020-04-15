//! Solves [Water Tree](http://codeforces.com/contest/343/problem/D).
//! To make a self-contained file for contest submission, dump each desired
//! module's contents directly here instead of the use statements.
//! Also, use the commented code in main() to employ standard I/O.
extern crate contest_algorithms;
use contest_algorithms::graph::Graph;
use contest_algorithms::range_query::{specs::AssignSum, StaticArq};
use contest_algorithms::scanner::Scanner;
use std::io;

const SAMPLE_INPUT: &[u8] = b"\
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
const SAMPLE_OUTPUT: &[u8] = b"\
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

fn solve<R: io::BufRead, W: io::Write>(scan: &mut Scanner<R>, out: &mut W) {
    let n = scan.token::<usize>();
    let mut tree = Graph::new(n, 2 * (n - 1));
    for _ in 1..n {
        let u = scan.token::<usize>() - 1;
        let v = scan.token::<usize>() - 1;
        tree.add_undirected_edge(u, v);
    }

    let mut l = vec![0; n];
    let mut r = vec![0; n];
    let mut p = vec![0; n];
    dfs(&tree, 0, &mut l, &mut r, &mut p, &mut 0);

    let mut arq = StaticArq::<AssignSum>::new(&vec![0; n + 1]);
    let q = scan.token::<usize>();
    for _ in 0..q {
        let c = scan.token::<usize>();
        let v = scan.token::<usize>() - 1;
        let len = (r[v] - l[v] + 1) as i64;
        let sum = arq.query(l[v], r[v]);
        if c == 1 {
            if sum != len {
                arq.update(p[v], p[v], &0);
                arq.update(l[v], r[v], &1);
            }
        } else if c == 2 {
            arq.update(l[v], l[v], &0);
        } else {
            let ans = if sum == len { 1 } else { 0 };
            writeln!(out, "{}", ans).ok();
        }
    }
}

#[test]
fn main() {
    let mut scan = Scanner::new(SAMPLE_INPUT);
    let mut out = vec![];
    solve(&mut scan, &mut out);

    assert_eq!(out, SAMPLE_OUTPUT);
}
