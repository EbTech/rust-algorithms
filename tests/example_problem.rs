/*extern crate algorithms;

use algorithms::scanner::*;
use algorithms::graph::*;

#[test]
fn integration() {
    let mut scan = Scanner::new();
    let n = scan.next::<usize>();
    let mut tree = Graph::new(n, n-1);
    for e in 0..n-1 {
        let u = scan.next::<usize>() - 1;
        let v = scan.next::<usize>() - 1;
        tree.add_edge(u, v);
    }
}*/
