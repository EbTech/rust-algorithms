//! Graph connectivity structures.
use super::Graph;

/// Helper struct that carries data needed for the depth-first searches in
/// ConnectivityGraph's constructor.
struct ConnectivityData {
    time: usize,
    vis: Box<[usize]>,
    low: Box<[usize]>,
    v_stack: Vec<usize>,
    e_stack: Vec<usize>,
}

impl ConnectivityData {
    fn new(num_v: usize) -> Self {
        Self {
            time: 0,
            vis: vec![0; num_v].into_boxed_slice(),
            low: vec![0; num_v].into_boxed_slice(),
            v_stack: Vec::new(),
            e_stack: Vec::new(),
        }
    }

    fn visit(&mut self, u: usize) {
        self.time += 1;
        self.vis[u] = self.time;
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
/// - 2-edge-connected components (2ECC),
/// - 2-vertex-connected components (2VCC)
///
/// Multiple-edges and self-loops are correctly handled.
pub struct ConnectivityGraph<'a> {
    // Immutable graph, frozen for the lifetime of the ConnectivityGraph object.
    pub graph: &'a Graph,
    /// ID of a vertex's CC, SCC or 2ECC, whichever applies.
    pub cc: Vec<usize>,
    /// ID of an edge's 2VCC, where applicable.
    pub vcc: Vec<usize>,
    /// Total number of CCs, SCCs or 2ECCs, whichever applies.
    pub num_cc: usize,
    /// Total number of 2VCCs, where applicable.
    pub num_vcc: usize,
}

impl<'a> ConnectivityGraph<'a> {
    /// Computes CCs (connected components), SCCs (strongly connected
    /// components), 2ECCs (2-edge-connected components), and/or 2VCCs
    /// (2-vertex-connected components), depending on the parameter and graph:
    /// - is_directed == true on directed graph: SCCs in rev-topological order
    /// - is_directed == true on undirected graph: CCs
    /// - is_directed == false on undirected graph: 2ECCs and 2VCCs
    /// - is_directed == false on directed graph: undefined behavior
    pub fn new(graph: &'a Graph, is_directed: bool) -> Self {
        let mut connect = Self {
            graph,
            cc: vec![0; graph.num_v()],
            vcc: vec![0; graph.num_e()],
            num_cc: 0,
            num_vcc: 0,
        };
        let mut data = ConnectivityData::new(graph.num_v());
        for u in 0..graph.num_v() {
            if data.vis[u] == 0 {
                if is_directed {
                    connect.scc(&mut data, u);
                } else {
                    connect.bcc(&mut data, u, graph.num_e() + 1);
                }
            }
        }
        connect
    }

    fn scc(&mut self, data: &mut ConnectivityData, u: usize) {
        data.visit(u);
        for (_, v) in self.graph.adj_list(u) {
            if data.vis[v] == 0 {
                self.scc(data, v);
            }
            if self.cc[v] == 0 {
                data.lower(u, data.low[v]);
            }
        }
        if data.vis[u] == data.low[u] {
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
            }).collect()
    }

    /// Gets the vertices of a graph according to a topological order of the
    /// strongly connected components. Most often used on DAGs.
    pub fn topological_sort(&self) -> Vec<usize> {
        let mut vertices = (0..self.graph.num_v()).collect::<Vec<_>>();
        vertices.sort_unstable_by_key(|&u| self.num_cc - self.cc[u]);
        vertices
    }

    fn bcc(&mut self, data: &mut ConnectivityData, u: usize, par: usize) {
        data.visit(u);
        for (e, v) in self.graph.adj_list(u) {
            if data.vis[v] == 0 {
                data.e_stack.push(e);
                self.bcc(data, v, e);
                data.lower(u, data.low[v]);
                if data.vis[u] <= data.low[v] {
                    // u is a cut vertex unless it's a one-child root
                    self.num_vcc += 1;
                    while let Some(top_e) = data.e_stack.pop() {
                        self.vcc[top_e] = self.num_vcc;
                        self.vcc[top_e ^ 1] = self.num_vcc;
                        if e ^ top_e <= 1 {
                            break;
                        }
                    }
                }
            } else if data.vis[v] < data.vis[u] && e ^ par != 1 {
                data.lower(u, data.vis[v]);
                data.e_stack.push(e);
            } else if v == u {
                // e is a self-loop
                self.num_vcc += 1;
                self.vcc[e] = self.num_vcc;
                self.vcc[e ^ 1] = self.num_vcc;
            }
        }
        if data.vis[u] == data.low[u] {
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

    /// In an undirected graph, determines whether u is an articulation vertex.
    pub fn is_cut_vertex(&self, u: usize) -> bool {
        if let Some(first_e) = self.graph.first[u] {
            self.graph
                .adj_list(u)
                .any(|(e, _)| self.vcc[first_e] != self.vcc[e])
        } else {
            false
        }
    }

    /// In an undirected graph, determines whether e is a bridge
    pub fn is_cut_edge(&self, e: usize) -> bool {
        let u = self.graph.endp[e ^ 1];
        let v = self.graph.endp[e];
        self.cc[u] != self.cc[v]
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_toposort() {
        let mut graph = Graph::new(4, 4);
        graph.add_edge(0, 2);
        graph.add_edge(3, 2);
        graph.add_edge(3, 1);
        graph.add_edge(1, 0);

        assert_eq!(
            ConnectivityGraph::new(&graph, true).topological_sort(),
            vec![3, 1, 0, 2]
        );
    }

    #[test]
    fn test_two_sat() {
        let mut graph = Graph::new(6, 8);
        let (x, y, z) = (0, 2, 4);

        graph.add_two_sat_clause(x, z);
        graph.add_two_sat_clause(y ^ 1, z ^ 1);
        graph.add_two_sat_clause(y, y);
        assert_eq!(
            ConnectivityGraph::new(&graph, true).two_sat_assign(),
            Some(vec![true, true, false])
        );

        graph.add_two_sat_clause(z, z);
        assert_eq!(ConnectivityGraph::new(&graph, true).two_sat_assign(), None);
    }

    #[test]
    fn test_biconnected() {
        let mut graph = Graph::new(3, 6);
        graph.add_undirected_edge(0, 1);
        graph.add_undirected_edge(1, 2);
        graph.add_undirected_edge(1, 2);

        let cg = ConnectivityGraph::new(&graph, false);
        let bridges = (0..graph.num_e())
            .filter(|&e| cg.is_cut_edge(e))
            .collect::<Vec<_>>();
        let articulation_points = (0..graph.num_v())
            .filter(|&u| cg.is_cut_vertex(u))
            .collect::<Vec<_>>();

        assert_eq!(bridges, vec![0, 1]);
        assert_eq!(articulation_points, vec![1]);
    }
}
