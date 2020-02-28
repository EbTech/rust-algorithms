use super::{DisjointSets, Graph};
use crate::graph::AdjListIterator;
use std::cmp::Reverse;

impl Graph {
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
        edges.sort_unstable_by_key(|&e| weights[e]);

        let mut components = DisjointSets::new(self.num_v());
        edges
            .into_iter()
            .filter(|&e| components.merge(self.endp[2 * e], self.endp[2 * e + 1]))
            .collect()
    }

    // Single-source shortest paths on a directed graph with non-negative weights
    pub fn dijkstra(&self, weights: &[u64], u: usize) -> Vec<u64> {
        assert_eq!(self.num_e(), weights.len());
        let mut dist = vec![u64::max_value(); weights.len()];
        let mut heap = std::collections::BinaryHeap::new();

        dist[u] = 0;
        heap.push((Reverse(0), 0));
        while let Some((Reverse(dist_u), u)) = heap.pop() {
            if dist[u] == dist_u {
                for (e, v) in self.adj_list(u) {
                    let dist_v = dist_u + weights[e];
                    if dist[v] > dist_v {
                        dist[v] = dist_v;
                        heap.push((Reverse(dist_v), v));
                    }
                }
            }
        }
        dist
    }

    pub fn dfs(&self, u: usize) -> DfsIterator {
        let adj_iters = (0..self.num_v())
            .map(|u| self.adj_list(u))
            .collect::<Vec<_>>();

        DfsIterator {
            visited: vec![false; self.num_v()],
            stack: vec![u],
            adj_iters,
        }
    }
}

pub struct DfsIterator<'a> {
    visited: Vec<bool>,
    stack: Vec<usize>,
    adj_iters: Vec<AdjListIterator<'a>>,
}

impl<'a> Iterator for DfsIterator<'a> {
    type Item = usize;

    /// Returns next vertex in the DFS
    fn next(&mut self) -> Option<Self::Item> {
        // Sources:
        // https://www.geeksforgeeks.org/iterative-depth-first-traversal/
        // https://en.wikipedia.org/wiki/Depth-first_search
        while let Some(&u) = self.stack.last() {
            if let Some((_, v)) = self.adj_iters[u].next() {
                if !self.visited[v] {
                    self.stack.push(v);
                }
            } else {
                self.stack.pop();
            }

            // Stack may contain same vertex twice. So
            // we return the popped item only
            // if it is not visited.
            if !self.visited[u] {
                self.visited[u] = true;
                return Some(u);
            }
        }

        None
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
        let mut graph = Graph::new(3, 6);
        graph.add_undirected_edge(0, 1);
        graph.add_undirected_edge(1, 2);
        graph.add_undirected_edge(2, 0);
        let weights = [7, 3, 5];

        let mst = graph.min_spanning_tree(&weights);
        let mst_cost = mst.iter().map(|&e| weights[e]).sum::<i64>();
        assert_eq!(mst, vec![1, 2]);
        assert_eq!(mst_cost, 8);
    }

    #[test]
    fn test_dijkstra() {
        let mut graph = Graph::new(3, 3);
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(2, 0);
        let weights = [7, 3, 5];

        let dist = graph.dijkstra(&weights, 0);
        assert_eq!(dist, vec![0, 7, 10]);
    }

    #[test]
    fn test_dfs() {
        let mut graph = Graph::new(4, 8);
        graph.add_edge(0, 2);
        graph.add_edge(2, 0);
        graph.add_edge(1, 2);
        graph.add_edge(0, 1);
        graph.add_edge(3, 3);
        graph.add_edge(2, 3);

        let dfs_search = graph.dfs(2).collect::<Vec<_>>();
        assert_eq!(dfs_search, vec![2, 3, 0, 1]);
    }

    #[test]
    fn test_dfs2() {
        let mut graph = Graph::new(5, 8);
        graph.add_edge(0, 2);
        graph.add_edge(2, 1);
        graph.add_edge(1, 0);
        graph.add_edge(0, 3);
        graph.add_edge(3, 4);
        graph.add_edge(4, 0);

        let dfs_search = graph.dfs(0).collect::<Vec<_>>();
        //Note this is not the only valid DFS
        assert_eq!(dfs_search, vec![0, 3, 4, 2, 1]);
    }

    #[test]
    fn test_dfs_space_complexity() {
        let num_v = 20;
        let mut graph = Graph::new(num_v, 0);
        for i in 0..num_v {
            for j in 0..num_v {
                graph.add_undirected_edge(i, j);
            }
        }

        let mut dfs_search = graph.dfs(7);
        let mut dfs_check = vec![];
        for _ in 0..num_v {
            dfs_check.push(dfs_search.next().unwrap());
            assert!(dfs_search.stack.len() <= num_v + 1);
        }

        dfs_check.sort();
        dfs_check.dedup();
        assert_eq!(0, dfs_check[0]);
        assert_eq!(num_v, dfs_check.len());
        assert_eq!(num_v - 1, dfs_check[num_v - 1]);
    }
}
