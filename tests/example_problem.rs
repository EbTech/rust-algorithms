// To make a self-contained file for contest submission, dump each desired
// module's contents directly here instead of these use statements.
extern crate algorithms;
use algorithms::scanner::*;
use algorithms::graph::*;

fn main1() {
    let stdin = std::io::stdin();
    let mut scan = Scanner::new(stdin.lock());
    let mut graph = Graph::new(1, 0);
    for _ in 0..0 {
        let u = scan.next::<usize>();
        let v = scan.next::<usize>();
        graph.add_edge(u, v);
    }
}

#[test]
fn main() {
    // If your contest solution requires a lot of stack space, make sure to
    // run it in a custom thread.
    std::thread::Builder::new()
        .stack_size(50_000_000)
        .spawn(main1)
        .unwrap()
        .join()
        .unwrap();
}
