//! Graph connectivity structures.

use super::graph::{DirectedGraph, UndirectedGraph};

/// Helper struct that carries data needed for the depth-first searches in
/// ConnectivityGraph's constructor.
struct ConnectivityData {
    time: usize,
    visited: Box<[usize]>,
    low: Box<[usize]>,
    v_stack: Vec<usize>,
    e_stack: Vec<usize>,
}

impl ConnectivityData {
    fn new(num_v: usize) -> Self {
        Self {
            time: 0,
            visited: vec![0; num_v].into_boxed_slice(),
            low: vec![0; num_v].into_boxed_slice(),
            v_stack: vec![],
            e_stack: vec![],
        }
    }

    fn visit(&mut self, u: usize) {
        self.time += 1;
        self.visited[u] = self.time;
        self.low[u] = self.time;
        self.v_stack.push(u);
    }

    fn lower(&mut self, u: usize, val: usize) {
        if self.low[u] > val {
            self.low[u] = val
        }
    }
}

/// Represents the decomposition of a graph into any of its constituent parts:
///
/// - Connected components (CC),
/// - Strongly connected components (SCC),
///
/// Multiple-edges and self-loops are correctly handled.
pub struct ConnectivityDirectedGraph<'a> {
    /// Immutable graph, frozen for the lifetime of the ConnectivityGraph object.
    pub graph: &'a DirectedGraph,
    /// ID of a vertex's CC, SCC or 2ECC, whichever applies. Range 1 to num_cc.
    pub cc: Vec<usize>,
    /// Total number of CCs, SCCs or 2ECCs, whichever applies.
    pub num_cc: usize,
}

impl<'a> ConnectivityDirectedGraph<'a> {
    /// Computes CCs (connected components), SCCs (strongly connected
    /// components),depending on the parameter and graph:
    pub fn new(graph: &'a DirectedGraph) -> Self {
        let mut connect = Self {
            graph,
            cc: vec![0; graph.num_v()],
            num_cc: 0,
        };
        let mut data = ConnectivityData::new(graph.num_v());
        for u in 0..graph.num_v() {
            if data.visited[u] == 0 {
                connect.scc(&mut data, u);
            }
        }
        connect
    }

    fn scc(&mut self, data: &mut ConnectivityData, u: usize) {
        data.visit(u);
        for (_, v) in self.graph.adj_list(u) {
            if data.visited[*v] == 0 {
                self.scc(data, *v);
            }
            if self.cc[*v] == 0 {
                data.lower(u, data.low[*v]);
            }
        }
        if data.visited[u] == data.low[u] {
            self.num_cc += 1;
            while let Some(v) = data.v_stack.pop() {
                self.cc[v] = self.num_cc;
                if v == u {
                    break;
                }
            }
        }
    }

    /// From the directed implication graph corresponding to a 2-SAT clause,
    /// finds a satisfying assignment if it exists or returns None otherwise.
    /// only for directed graphs
    pub fn two_sat_assign(&self) -> Option<Vec<bool>> {
        (0..self.graph.num_v() / 2)
            .map(|i| {
                let scc_true = self.cc[2 * i];
                let scc_false = self.cc[2 * i + 1];
                if scc_true == scc_false {
                    None
                } else {
                    Some(scc_true < scc_false)
                }
            })
            .collect()
    }

    /// Gets the vertices of a graph according to a topological order of the
    /// strongly connected components. Most often used on DAGs.
    pub fn topological_sort(&self) -> Vec<usize> {
        let mut vertices = (0..self.graph.num_v()).collect::<Vec<_>>(); //;turbofish
        vertices.sort_unstable_by_key(|&u| self.num_cc - self.cc[u]);
        vertices
    }
}

/// Represents the decomposition of a graph into any of its constituent parts:
///
/// - Connected components (CC),
/// - 2-edge-connected components (2ECC),
/// - 2-vertex-connected components (2VCC)
///
/// Multiple-edges and self-loops are correctly handled.
/// this class find connected components, as well as 2 connected components
pub struct ConnectivityUndirectedGraph<'a> {
    /// Immutable undiredted graph, frozen for the lifetime of this object.
    pub graph: &'a UndirectedGraph,
    /// ID of a vertex's CC or 2ECC, whichever applies. Range 1 to num_cc.
    pub cc: Vec<usize>,
    /// ID of an edge's 2VCC, where applicable. Ranges from 1 to num_vcc.
    pub vcc: Vec<usize>,
    /// keeps track of articuallation points, aka cut vertex's
    pub is_articulation_point: Vec<bool>,
    /// Total number of CCs, SCCs or 2ECCs, whichever applies.
    pub num_cc: usize,
    /// Total number of 2VCCs, where applicable.
    pub num_vcc: usize,
}

impl<'a> ConnectivityUndirectedGraph<'a> {
    /// Computes CCs (connected components), SCCs (strongly connected
    /// components), 2ECCs (2-edge-connected components), and/or 2VCCs
    /// (2-vertex-connected components), depending on the parameter and graph:
    /// - is_directed == true on directed graph: SCCs in rev-topological order
    /// - is_directed == true on undirected graph: CCs
    /// - is_directed == false on undirected graph: 2ECCs and 2VCCs
    /// - is_directed == false on directed graph: undefined behavior
    pub fn new(graph: &'a UndirectedGraph) -> Self {
        let mut connect = Self {
            graph,
            cc: vec![0; graph.num_v()],
            vcc: vec![0; graph.num_e()],
            is_articulation_point: vec![false; graph.num_v()],
            num_cc: 0,
            num_vcc: 0,
        };
        let mut data = ConnectivityData::new(graph.num_v());
        for u in 0..graph.num_v() {
            if data.visited[u] == 0 {
                connect.biconnected(&mut data, u, usize::MAX);
            }
        }
        connect
    }

    ///Tarjans algorithm  for finding cut vertex. this also find biconnected components
    /// https://www.geeksforgeeks.org/articulation-points-or-cut-vertices-in-a-graph/
    fn biconnected(&mut self, data: &mut ConnectivityData, u: usize, parent: usize) {
        data.visit(u);
        let mut children = 0usize;
        for (er, vr) in self.graph.adj_list(u) {
            let e = *er;
            let v = *vr;
            if data.visited[v] == 0 {
                children += 1;
                data.e_stack.push(e);
                self.biconnected(data, v, u);

                data.lower(u, data.low[v]);
                if parent < usize::MAX && data.visited[u] <= data.low[v] {
                    // u is a cut vertex unless it's a one-child root
                    self.num_vcc += 1;
                    //data.e_stack.pop();
                    while let Some(top_e) = data.e_stack.pop() {
                        if self.vcc[top_e] == 0 {
                            //never been assigned b4
                            self.vcc[top_e] = self.num_vcc;
                        }
                        if e == top_e {
                            break;
                        }
                    }
                    self.is_articulation_point[u] = true;
                }
            } else if data.visited[v] < data.visited[u] &&  e != parent{
                data.lower(u, data.visited[v]);
                data.e_stack.push(e);
            } else if v == u {
                // e is a self-loop
                self.num_vcc += 1;
                self.vcc[e] = self.num_vcc;
            }
        }
        // if u is a root of dfs tree and has two or more children
        if parent == usize::MAX && children > 1 {
            self.num_vcc += 1;

            while let Some(top_e) = data.e_stack.pop() {
                self.vcc[top_e] = self.num_vcc;
            }

            self.is_articulation_point[u] = true;
        }
        if data.visited[u] == data.low[u] {
            // par is a cut edge unless par==-1
            self.num_cc += 1;
            while let Some(v) = data.v_stack.pop() {
                self.cc[v] = self.num_cc;
                if v == u {
                    break;
                }
            }
        }
    }

    /// In an undirected graph, determines whether u is an cut vertex.
    pub fn is_cut_vertex(&self, u: usize) -> bool {
        //return self.is_articulation_point[u];

        let first_e = self.graph.adj_lists[u][0].0;
        self.graph
            .adj_list(u)
            .any(|(e, _)| self.vcc[first_e] != self.vcc[*e])
    }

    /// In an undirected graph, determines whether e is a bridge
    pub fn is_cut_edge(&self, e: usize) -> bool {
        let ev = Vec::from_iter(&self.graph.edges[e]);
        self.cc[*ev[0]] != self.cc[*ev[1]]
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_toposort() {
        let mut graph = DirectedGraph::new(4, 5);
        graph.add_edge(0, 0);
        graph.add_edge(0, 2);
        graph.add_edge(3, 2);
        graph.add_edge(3, 1);
        graph.add_edge(1, 0);

        assert_eq!(
            ConnectivityDirectedGraph::new(&graph).topological_sort(),
            vec![3, 1, 0, 2]
        );
    }

    #[test]
    fn test_two_sat() {
        let mut graph = DirectedGraph::new(6, 8);
        let (x, y, z) = (0, 2, 4);

        graph.add_two_sat_clause(x, z);
        graph.add_two_sat_clause(y ^ 1, z ^ 1);
        graph.add_two_sat_clause(y, y);
        assert_eq!(
            ConnectivityDirectedGraph::new(&graph).two_sat_assign(),
            Some(vec![true, true, false])
        );

        graph.add_two_sat_clause(z, z);
        assert_eq!(
            ConnectivityDirectedGraph::new(&graph).two_sat_assign(),
            None
        );
    }

    #[test]
    fn test_biconnected() {
        let mut graph = UndirectedGraph::new(3, 6);
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(1, 2);

        let cg = ConnectivityUndirectedGraph::new(&graph);
        let bridges = (0..graph.num_e())
            .filter(|&e| cg.is_cut_edge(e))
            .collect::<Vec<_>>();
        let articulation_points = (0..graph.num_v())
            .filter(|&u| cg.is_cut_vertex(u))
            .collect::<Vec<_>>();

        for idx in 0..graph.num_e(){

            if cg.is_cut_edge(idx){
                println!(" edge {} is a cut edge",idx)
            }
            else{
                println!(" edge {} is not a cut edge",idx)
            }
        }
        assert_eq!(bridges, vec![0]);
        assert_eq!(articulation_points, vec![1]);
    }

    #[test]
    fn test_articulation_points() {
        let mut graph = UndirectedGraph::new(7, 8);
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(0, 2);
        graph.add_edge(1, 3);
        graph.add_edge(3, 5);
        graph.add_edge(5, 4);
        graph.add_edge(4, 1);
        graph.add_edge(1, 6);

        let cg = ConnectivityUndirectedGraph::new(&graph);
       
        let articulation_points = (0..graph.num_v())
            .filter(|&u| cg.is_cut_vertex(u))
            .collect::<Vec<_>>();

        assert_eq!(articulation_points, vec![1]);
    }
}
