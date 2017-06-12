use graph::Graph;
use ::std::cmp::min;

// Represents the decomposition of a graph into any of its:
// - Connected components (CC),
// - Strongly connected components (SCC),
// - 2-edge-connected components (2ECC),
// - 2-vertex-connected components (2VCC)
// Multiple-edges and self-loops should be correctly handled.
pub struct ConnectivityGraph<'a> {
    pub graph: &'a Graph,
    pub cc: Vec<usize>,  // stores id of a vertex's CC, SCC or 2ECC
    pub vcc: Vec<usize>, // stores id of an edge's 2VCC
    pub num_cc: usize,
    pub num_vcc: usize
}

impl<'a> ConnectivityGraph<'a> {
    // Computes the SCCs of a directed graph in reverse topological order, or
    // the 2ECCs/2VCCS of an undirected graph. Can also get CCs by passing an
    // undirected graph with is_directed == true.
    pub fn new(graph: &'a Graph, is_directed: bool) -> ConnectivityGraph {
        let mut connect = ConnectivityGraph {
            graph: graph,
            cc: vec![0; graph.num_v()],
            vcc: vec![0; graph.num_e()],
            num_cc: 0,
            num_vcc: 0
        };
        let mut t = 0;
        let mut vis = vec![0; graph.num_v()];
        let mut low = vec![0; graph.num_v()];
        let mut verts = Vec::new();
        let mut edges = Vec::new();
        for u in 0..graph.num_v() {
            if vis[u] == 0 {
                if is_directed {
                    connect.scc(u, &mut t, &mut vis, &mut low, &mut verts);
                }
                else {
                    connect.bcc(u, graph.num_e() + 1, &mut t, &mut vis,
                                &mut low, &mut verts, &mut edges);
                }
            }
        }
        connect
    }
    
    fn scc(&mut self, u: usize, t: &mut usize, vis: &mut [usize],
           low: &mut [usize], verts: &mut Vec<usize>) {
        *t += 1;
        vis[u] = *t;
        low[u] = *t;
        verts.push(u);
        for (_, v) in self.graph.adj_list(u) {
            if vis[v] == 0 { self.scc(v, t, vis, low, verts); }
            if self.cc[v] == 0 { low[u] = min(low[u], low[v]); }
        }
        if vis[u] == low[u] {
            self.num_cc += 1;
            while let Some(v) = verts.pop() {
                self.cc[v] = self.num_cc;
                if v == u { break; }
            }
        }
    }
    
    // From the directed implication graph corresponding to a 2-SAT clause,
    // finds a satisfying assignment if it exists or returns None otherwise.
    pub fn two_sat_assign(&self) -> Option<Vec<bool>> {
        (0..self.graph.num_v()/2).map( |i| {
            let scc_true = self.cc[2*i];
            let scc_false = self.cc[2*i+1];
            if scc_true == scc_false { None } else { Some(scc_true < scc_false) }
        }).collect()
    }
    
    // Gets the vertices of a directed acyclic graph (DAG) in topological order.
    pub fn topological_sort(&self) -> Vec<usize> {
        let mut vertices = (0..self.graph.num_v()).collect::<Vec<_>>();
        vertices.sort_by_key(|&u| self.num_cc - self.cc[u]);
        vertices
    }
    
    fn bcc(&mut self, u: usize, par: usize, t: &mut usize, vis: &mut [usize],
           low: &mut [usize], verts: &mut Vec<usize>, edges: &mut Vec<usize>) {
        *t += 1;
        vis[u] = *t;
        low[u] = *t;
        verts.push(u);
        for (e, v) in self.graph.adj_list(u) {
          if vis[v] == 0 {
              edges.push(e);
              self.bcc(v, e, t, vis, low, verts, edges);
              low[u] = min(low[u], low[v]);
              if vis[u] <= low[v] { // u is a cut vertex unless it's a one-child root
                  self.num_vcc += 1;
                  while let Some(top_e) = edges.pop() {
                      self.vcc[top_e] = self.num_vcc;
                      self.vcc[top_e ^ 1] = self.num_vcc;
                      if e ^ top_e <= 1 { break; }
                  }
              }
          }
          else if vis[v] < vis[u] && e ^ par != 1 {
              low[u] = min(low[u], vis[v]);
              edges.push(e);
          }
          else if v == u { // e is a self-loop
              self.num_vcc += 1;
              self.vcc[e] = self.num_vcc;
              self.vcc[e ^ 1] = self.num_vcc;
          }
        } 
        if vis[u] == low[u] { // par is a cut edge unless par==-1
            self.num_cc += 1;
            while let Some(v) = verts.pop() {
                self.cc[v] = self.num_cc;
                if v == u { break; }
            } 
        }
    }
    
    // In an undirected graph, determines whether u is an articulation vertex.
    pub fn is_cut_vertex(&self, u: usize) -> bool {
        if let Some(first_e) = self.graph.first[u] {
            self.graph.adj_list(u).any(|(e, _)| {self.vcc[first_e] != self.vcc[e]})
        }
        else {
            false
        }
    }
    
    // In an undirected graph, determines whether v is a bridge
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
    fn test_toposort()
    {
        let mut graph = Graph::new(4, 4);
        graph.add_edge(0, 2);
        graph.add_edge(3, 2);
        graph.add_edge(3, 1);
        graph.add_edge(1, 0);
        
        assert_eq!(ConnectivityGraph::new(&graph, true).topological_sort(),
                   vec![3, 1, 0, 2]);
    }
    
    #[test]
    fn test_two_sat()
    {
        let mut graph = Graph::new(6, 8);
        let (x, y, z) = (0, 2, 4);
        
        graph.add_two_sat_clause(x, z);
        graph.add_two_sat_clause(y^1, z^1);
        graph.add_two_sat_clause(y, y);
        assert_eq!(ConnectivityGraph::new(&graph, true).two_sat_assign(),
                   Some(vec![true, true, false]));
            
        graph.add_two_sat_clause(z, z);
        assert_eq!(ConnectivityGraph::new(&graph, true).two_sat_assign(), None);
    }
    
    #[test]
    fn test_biconnected()
    {
        let mut graph = Graph::new(3, 6);
        graph.add_undirected_edge(0, 1);
        graph.add_undirected_edge(1, 2);
        graph.add_undirected_edge(1, 2);
        
        let cg = ConnectivityGraph::new(&graph, false);
        let bridges = (0..graph.num_e())
                    .filter(|&e| cg.is_cut_edge(e)).collect::<Vec<_>>();
        let articulation_points = (0..graph.num_v())
                    .filter(|&u| cg.is_cut_vertex(u)).collect::<Vec<_>>();
        
        assert_eq!(bridges, vec![0, 1]);
        assert_eq!(articulation_points, vec![1]);
    }
}
