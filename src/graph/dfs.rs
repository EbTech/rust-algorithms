use super::Graph;
use crate::graph::AdjListIterator;
use bit_vec::BitVec;

impl Graph {
    pub fn dfs(&self, v: usize) -> DfsIterator {
        // Create a stack for DFS
        let mut stack: Vec<usize> = Vec::new();

        let adj_iters = (0..self.num_v())
            .map(|u| self.adj_list(u))
            .collect::<Vec<_>>();

        // Push the current source node.
        stack.push(v);

        DfsIterator {
            visited: BitVec::from_elem(self.num_v(), false),
            stack,
            adj_iters,
        }
    }
}
pub struct DfsIterator<'a> {
    //is vertex visited
    visited: BitVec,
    //stack of vertices
    stack: Vec<usize>,
    adj_iters: Vec<AdjListIterator<'a>>,
}

impl<'a> Iterator for DfsIterator<'a> {
    type Item = usize;

    /// Returns next vertex in the DFS
    fn next(&mut self) -> Option<Self::Item> {
        //Sources:
        // https://www.geeksforgeeks.org/iterative-depth-first-traversal/
        // https://en.wikipedia.org/wiki/Depth-first_search
        while let Some(&s) = self.stack.last() {
            //Does s still have neighbors we need to process?
            if let Some((_, s_nbr)) = self.adj_iters[s].next() {
                if !self.visited[s_nbr] {
                    self.stack.push(s_nbr);
                }
            } else {
                //s has no more neighbors, we can pop it off the stack
                self.stack.pop();
            }

            // Stack may contain same vertex twice. So
            // we return the popped item only
            // if it is not visited.
            if !self.visited[s] {
                self.visited.set(s, true);
                return Some(s);
            }
        }

        None
    }
}

#[cfg(test)]
mod test {
    use super::*;

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
