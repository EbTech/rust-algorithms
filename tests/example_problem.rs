// To make a single-file Codeforces contest submission, dump the module contents
// directly here instead of these use statements.
extern crate algorithms;
use algorithms::scanner::*;
use algorithms::graph::*;

fn main1() {
    let mut scan = Scanner::new();
    let mut graph = Graph::new(1, 0);
    for _ in 0..0 {
        let u = scan.next::<usize>();
        let v = scan.next::<usize>();
        graph.add_edge(u, v);
    }
}

#[test]
fn integration() {
    // If your contest solution requires a lot of stack space, make sure to
    // run it on a custom thread.
    std::thread::Builder::new().stack_size(50_000_000)
        .spawn(main1).unwrap().join().unwrap();
}
