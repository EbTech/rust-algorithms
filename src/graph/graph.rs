//! Basic graph module without explicit support for deletion.
//!
//! # Panics
//!
//! All methods will panic if given an out-of-bounds element index.
/// A compact graph representation. Edges are numbered in order of insertion.
/// Each adjacency list consists of all edges pointing out from a given vertex.
///
use std::collections::HashMap;

pub struct DirectedGraph {
    /// Maps a vertex id to the first edge in its adjacency list.
    pub first: Vec<Option<usize>>,
    /// Maps an edge id to the next edge in the same adjacency list.
    pub next: Vec<Option<usize>>,
    /// Maps an edge id to the vertex that it points to.
    pub endp: Vec<usize>,
    /// Set containing all the edges, used for quick look up
    pub edge_weights: HashMap<(usize, usize), i64>,
}

impl DirectedGraph {
    /// Initializes a graph with vmax vertices and no edges. To reduce
    /// unnecessary allocations, emax_hint should be close to the number of
    /// edges that will be inserted.
    pub fn new(vmax: usize, emax_hint: usize) -> Self {
        Self {
            first: vec![None; vmax],
            next: Vec::with_capacity(emax_hint),
            endp: Vec::with_capacity(emax_hint),
            edge_weights: HashMap::new(),
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
        self.edge_weights.insert((u, v), 1i64);
    }

    /// Adds a weighted directed edge from u to v.
    pub fn add_weighted_edge(&mut self, u: usize, v: usize, w: i64) {
        self.next.push(self.first[u]);
        self.first[u] = Some(self.num_e());
        self.endp.push(v);
        self.edge_weights.insert((u, v), w);
    }

    /// this retrieves a weight vector, where index is the edge index
    /// probably this should not be public, since edge index is internal representation
    pub fn get_weights(&self) -> Vec<i64> {
        let mut ret = vec![0i64; self.num_e()];

        for idx in 0..self.num_v() {
            for e_iter in self.adj_list(idx) {
                ret[e_iter.0] = self.edge_weights[&(idx, e_iter.1)];
            }
        }

        ret
    }

    /// If we think of each even-numbered vertex as a variable, and its
    /// odd-numbered successor as its negation, then we can build the
    /// implication graph corresponding to any 2-CNF formula.
    /// Note that u||v == !u -> v == !v -> u.
    pub fn add_two_sat_clause(&mut self, u: usize, v: usize) {
        self.add_edge(u ^ 1, v);
        self.add_edge(v ^ 1, u);
    }

    /// This tests if an edge is contained here
    pub fn has_edge(&self, u: usize, v: usize) -> bool {
        self.edge_weights.contains_key(&(u, v))
    }

    /// Gets vertex u's adjacency list.
    pub fn adj_list(&self, u: usize) -> AdjListIterator {
        AdjListIterator {
            graph: self,
            next_e: self.first[u],
        }
    }
}

/// An iterator for convenient adjacency list traversal.
pub struct AdjListIterator<'a> {
    graph: &'a DirectedGraph,
    next_e: Option<usize>,
}

impl<'a> Iterator for AdjListIterator<'a> {
    /// first is vertex index, second is edge index
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

pub struct UndirectedGraph {
    /// underlying representation. wouldnt it be nice if we had inheritance
    pub directed_graph: DirectedGraph,
}

impl UndirectedGraph {
    /// Initializes a graph with vmax vertices and no edges. To reduce
    /// unnecessary allocations, emax_hint should be close to the number of
    /// edges that will be inserted.
    pub fn new(vmax: usize, emax_hint: usize) -> Self {
        Self {
            directed_graph: DirectedGraph::new(vmax, 2*emax_hint),
        }
    }

    /// Returns the number of vertices.
    pub fn num_v(&self) -> usize {
        self.directed_graph.num_v()
    }

    /// Returns the number of edges, double-counting Undirected edges.
    pub fn num_e(&self) -> usize {
        self.directed_graph.num_e()
    }

    /// Adds a directed edge from u to v.
    pub fn add_edge(&mut self, u: usize, v: usize) {
        self.directed_graph.add_edge(u, v);
        self.directed_graph.add_edge(v, u);
    }

    /// Adds a weighted directed edge from u to v.
    pub fn add_weighted_edge(&mut self, u: usize, v: usize, w: i64) {
        self.directed_graph.add_weighted_edge(u, v, w);
        self.directed_graph.add_weighted_edge(v, u, w);
    }

    /// this retrieves a weight vector, where index is the edge index
    /// probably this should not be public, since edge index is internal representation
    pub fn get_weights(&self) -> Vec<i64> {
        self.directed_graph
            .get_weights()
            .into_iter()
            .enumerate()
            .filter(|e_idx| e_idx.0 % 2 == 0)
            .map(|(_, v)| v)
            .collect::<Vec<_>>() //turbofish!
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_adj_list() {
        let mut graph = DirectedGraph::new(5, 6);
        graph.add_edge(2, 3);
        graph.add_edge(2, 4);
        graph.add_edge(4, 1);
        graph.add_edge(1, 2);
        graph.add_edge(0, 2);
        graph.add_edge(2, 0);

        let adj = graph.adj_list(2).collect::<Vec<_>>();

        assert_eq!(adj, vec![(5, 0), (1, 4), (0, 3)]);
        for (e, v) in adj {
            assert_eq!(v, graph.endp[e]);
        }
    }
}
