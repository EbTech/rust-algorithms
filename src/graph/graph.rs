//! Basic graph module without explicit support for deletion.
//!
//! # Panics
//!
//! All methods will panic if given an out-of-bounds element index.
use core::slice::Iter;
/// A compact graph representation. Edges are numbered in order of insertion.
/// Each adjacency list consists of all edges pointing out from a given vertex.
///
use std::collections::HashMap;

pub type AdjListIterator<'a> = Iter<'a, (usize, usize)>;

pub struct DirectedGraph {
    /// adjacency list. each vertex has a list of (edge index, destination vertex index)
    pub adj_lists: Vec<Vec<(usize, usize)>>,
    /// Maps an edge id to the vertex that it points to.
    pub edges: Vec<(usize,usize)>,
    /// edge weights
    pub edge_weights: Vec<i64>,
}

impl DirectedGraph {
    /// Initializes a graph with vmax vertices and no edges. To reduce
    /// unnecessary allocations, emax_hint should be close to the number of
    /// edges that will be inserted.
    pub fn new(vmax: usize, emax_hint: usize) -> Self {
        Self {
            adj_lists: vec![Vec::with_capacity(vmax); vmax],
            edges: Vec::with_capacity(emax_hint),
            edge_weights: Vec::with_capacity(emax_hint),
        }
    }

    /// Returns the number of vertices.
    pub fn num_v(&self) -> usize {
        self.adj_lists.len()
    }

    /// Returns the number of edges, double-counting undirected edges.
    pub fn num_e(&self) -> usize {
        self.edges.len()
    }

    /// Adds a directed edge from u to v.
    pub fn add_edge(&mut self, u: usize, v: usize) {
        self.add_weighted_edge(u, v, 1i64);
    }

    /// Adds a weighted directed edge from u to v.
    pub fn add_weighted_edge(&mut self, u: usize, v: usize, w: i64) {
        self.edges.push((u,v));
        self.edge_weights.push(w);
        self.adj_lists[u].push((self.edges.len() - 1, v));
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
        self.adj_lists[u].iter()
    }
}

pub struct UndirectedGraph {
    /// adjacency list. each vertex has a list of (edge index, neighor vertex index)
    pub adj_lists: Vec<Vec<(usize, usize)>>,
    /// Maps an edge id to vertices. is stored as smalles index first
    pub edges: Vec<(usize, usize)>,
    /// edge weights
    pub edge_weights: Vec<i64>,
}

impl UndirectedGraph {
    /// Initializes a graph with vmax vertices and no edges. To reduce
    /// unnecessary allocations, emax_hint should be close to the number of
    /// edges that will be inserted.
    pub fn new(vmax: usize, emax_hint: usize) -> Self {
        Self {
            adj_lists: vec![Vec::with_capacity(vmax); vmax],
            edges: Vec::with_capacity(emax_hint),
            edge_weights: Vec::with_capacity(emax_hint),
        }
    }

    /// Returns the number of vertices.
    pub fn num_v(&self) -> usize {
        self.adj_lists.len()
    }

    /// Returns the number of edges, double-counting Undirected edges.
    pub fn num_e(&self) -> usize {
        self.edges.len()
    }

    /// Adds a directed edge from u to v.
    pub fn add_edge(&mut self, u: usize, v: usize) {
        self.add_weighted_edge(u, v, 1i64);
    }

    /// Adds a weighted directed edge from u to v.
    pub fn add_weighted_edge(&mut self, u: usize, v: usize, w: i64) {
        let minv = std::cmp::min(u, v);
        let maxv = std::cmp::max(u, v);
        self.edges.push((minv, maxv));
        self.edge_weights.push(w);
        self.adj_lists[u].push((self.edges.len() - 1, v));
        self.adj_lists[v].push((self.edges.len() - 1, u));
    }

    /// Gets vertex u's adjacency list.
    pub fn adj_list(&self, u: usize) -> AdjListIterator {
        self.adj_lists[u].iter()
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

        for (e, v) in adj {
            assert_eq!(*v, graph.edges[*e].1);
        }
    }

    #[test]
    fn test_undirected_graph_basic() {
        let mut graph = UndirectedGraph::new(5, 6);
        graph.add_edge(2, 3);
        graph.add_edge(4, 3);
        graph.add_edge(1, 2);
        graph.add_edge(1, 0);

        assert_eq!(4, graph.num_e());
        assert_eq!(5, graph.num_v());
    }
}
