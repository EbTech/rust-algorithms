//! Basic graph module without explicit support for deletion.
//!
//! # Panics
//!
//! All methods will panic if given an out-of-bounds element index.
pub mod flow;
pub mod connectivity;

/// Represents a union of disjoint sets. Each set's elements are arranged in a
/// tree, whose root is the set's representative.
pub struct DisjointSets {
    parent: Vec<usize>,
}

impl DisjointSets {
    /// Initializes disjoint sets containing one element each.
    pub fn new(size: usize) -> Self {
        Self { parent: (0..size).collect() }
    }

    /// Finds the set's representative. Do path compression along the way to make
    /// future queries faster.
    pub fn find(&mut self, u: usize) -> usize {
        let pu = self.parent[u];
        if pu != u {
            self.parent[u] = self.find(pu);
        }
        self.parent[u]
    }

    /// Merges the sets containing u and v into a single set containing their
    /// union. Returns true if u and v were previously in different sets.
    pub fn merge(&mut self, u: usize, v: usize) -> bool {
        let (pu, pv) = (self.find(u), self.find(v));
        self.parent[pu] = pv;
        pu != pv
    }
}

/// A compact graph representation. Edges are numbered in order of insertion.
pub struct Graph {
    first: Vec<Option<usize>>,
    next: Vec<Option<usize>>,
    endp: Vec<usize>,
}

impl Graph {
    /// Initializes a graph with vmax vertices and no edges. To reduce
    /// unnecessary allocations, emax_hint should be close to the number of
    /// edges that will be inserted.
    pub fn new(vmax: usize, emax_hint: usize) -> Self {
        Self {
            first: vec![None; vmax],
            next: Vec::with_capacity(emax_hint),
            endp: Vec::with_capacity(emax_hint),
        }
    }

    /// Returns the number of vertices.
    pub fn num_v(&self) -> usize {
        self.first.len()
    }

    /// Returns the number of edges, double-counting undirected edges.
    pub fn num_e(&self) -> usize {
        self.endp.len()
    }

    /// Adds a directed edge from u to v.
    pub fn add_edge(&mut self, u: usize, v: usize) {
        self.next.push(self.first[u]);
        self.first[u] = Some(self.num_e());
        self.endp.push(v);
    }

    /// An undirected edge is two directed edges. If edges are added only via
    /// this funcion, the reverse of any edge e can be found at e^1.
    pub fn add_undirected_edge(&mut self, u: usize, v: usize) {
        self.add_edge(u, v);
        self.add_edge(v, u);
    }

    /// If we think of each even-numbered vertex as a variable, and its
    /// odd-numbered successor as its negation, then we can build the
    /// implication graph corresponding to any 2-CNF formula.
    /// Note that u||v == !u -> v == !v -> u.
    pub fn add_two_sat_clause(&mut self, u: usize, v: usize) {
        self.add_edge(u ^ 1, v);
        self.add_edge(v ^ 1, u);
    }

    /// Gets vertex u's adjacency list.
    pub fn adj_list(&self, u: usize) -> AdjListIterator {
        AdjListIterator {
            graph: self,
            next_e: self.first[u],
        }
    }

    /// Finds the sequence of edges in an Euler path starting from u, assuming
    /// it exists and that the graph is directed. Undefined behavior if this
    /// precondition is violated. To extend this to undirected graphs, maintain
    /// a visited array to skip the reverse edge.
    pub fn euler_path(&self, u: usize) -> Vec<usize> {
        let mut adj_iters = (0..self.num_v())
            .map(|u| self.adj_list(u))
            .collect::<Vec<_>>();
        let mut edges = Vec::with_capacity(self.num_e());
        self.euler_recurse(u, &mut adj_iters, &mut edges);
        edges.reverse();
        edges
    }

    // Helper function used by euler_path. Note that we can't use a for-loop
    // that would consume the adjacency list as recursive calls may need it.
    fn euler_recurse(&self, u: usize, adj: &mut [AdjListIterator], edges: &mut Vec<usize>) {
        while let Some((e, v)) = adj[u].next() {
            self.euler_recurse(v, adj, edges);
            edges.push(e);
        }
    }

    /// Kruskal's minimum spanning tree algorithm on an undirected graph.
    pub fn min_spanning_tree(&self, weights: &[i64]) -> Vec<usize> {
        assert_eq!(self.num_e(), 2 * weights.len());
        let mut edges = (0..weights.len()).collect::<Vec<_>>();
        edges.sort_by_key(|&e| weights[e]);

        let mut components = DisjointSets::new(self.num_v());
        edges
            .into_iter()
            .filter(|&e| {
                components.merge(self.endp[2 * e], self.endp[2 * e + 1])
            })
            .collect()
    }
}

/// An iterator for convenient adjacency list traversal.
pub struct AdjListIterator<'a> {
    graph: &'a Graph,
    next_e: Option<usize>,
}

impl<'a> Iterator for AdjListIterator<'a> {
    type Item = (usize, usize);

    /// Produces an outgoing edge and vertex.
    fn next(&mut self) -> Option<Self::Item> {
        self.next_e.map(|e| {
            let v = self.graph.endp[e];
            self.next_e = self.graph.next[e];
            (e, v)
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_euler() {
        let mut graph = Graph::new(3, 4);
        graph.add_edge(0, 1);
        graph.add_edge(1, 0);
        graph.add_edge(1, 2);
        graph.add_edge(2, 1);

        assert_eq!(graph.euler_path(0), vec![0, 2, 3, 1]);
    }

    #[test]
    fn test_min_spanning_tree() {
        let mut graph = Graph::new(3, 3);
        graph.add_undirected_edge(0, 1);
        graph.add_undirected_edge(1, 2);
        graph.add_undirected_edge(2, 0);
        let weights = [7, 3, 5];

        let mst = graph.min_spanning_tree(&weights);
        let mst_cost = mst.iter().map(|&e| weights[e]).sum::<i64>();
        assert_eq!(mst, vec![1, 2]);
        assert_eq!(mst_cost, 8);
    }
}
