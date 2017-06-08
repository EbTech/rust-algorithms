

// Represents a union of disjoint sets. Each set's elements are arranged in a
// tree, whose root is the set's representative.
pub struct DisjointSets {
    parent: Vec<usize>
}

impl DisjointSets {
    // Initialize disjoint sets containing one element each.
    pub fn new(size: usize) -> DisjointSets {
        DisjointSets { parent: (0..size).collect() }
    }
    
    // Find the set's representative. Do path compression along the way to make
    // future queries faster.
    pub fn find(&mut self, u: usize) -> usize {
        let pu = self.parent[u];
        if pu != u { self.parent[u] = self.find(pu); }
        self.parent[u]
    }
    
    // Merge the sets containing u and v into a single set containing their
    // union. Returns true if u and v were previously in different sets.
    pub fn merge(&mut self, u: usize, v: usize) -> bool {
        let (pu, pv) = (self.find(u), self.find(v));
        self.parent[pu] = pv;
        pu != pv
    }
}

// A compact graph representation.
pub struct Graph {
    pub first: Vec<Option<usize>>,
    pub next: Vec<Option<usize>>,
    pub endp: Vec<usize>,
}

impl Graph {
    // Initialize a graph with vmax vertices and no edges. For best efficiency,
    // emax should be a tight upper bound on the number of edges to insert.
    pub fn new(vmax: usize, emax: usize) -> Graph {
        Graph {
            first: vec![None; vmax],
            next: Vec::with_capacity(emax),
            endp: Vec::with_capacity(emax)
        }
    }
    
    // Utilities to compute the number of vertices and edges.
    pub fn num_v(&self) -> usize { self.first.len() }
    pub fn num_e(&self) -> usize { self.next.len() }
    
    // Add a directed edge from u to v.
    pub fn add_edge(&mut self, u: usize, v: usize) {
        self.next.push(self.first[u]);
        self.first[u] = Some(self.endp.len());
        self.endp.push(v);
    }
    
    // An undirected edge is two directed edges. If edges are added only via
    // this funcion, the reverse of any edge e can be found at e^1.
    pub fn add_undirected_edge(&mut self, u: usize, v: usize) {
        self.add_edge(u, v);
        self.add_edge(v, u);
    }
    
    // If we think of each even-numbered vertex as a variable, and its successor
    // as its negation, then we can build the implication graph corresponding
    // to any 2-CNF formula. Note that u||v == !u -> v == !v -> u.
    pub fn add_two_sat_clause(&mut self, u: usize, v: usize) {
        self.add_edge(u^1, v);
        self.add_edge(v^1, u);
    }
    
    // Gets vertex u's adjacency list.
    pub fn adj_list<'a>(&'a self, u: usize) -> AdjListIterator<'a> {
        AdjListIterator {
            graph: self,
            next_e: self.first[u]
        }
    }
    
    // Helper function used by euler_path. Note that we can't consume the
    // adjacency list in a for loop because recursive calls may need it.
    fn euler_recurse(&self, u: usize, adj: &mut [AdjListIterator], edges: &mut Vec<usize>) {
        while let Some((e, v)) = adj[u].next() {
            self.euler_recurse(v, adj, edges);
            edges.push(e);
        }
    }
    
    // Finds the sequence of edges in an Euler path starting from u, assuming it
    // exists and that the graph is directed. To extend this to undirected
    // graphs, keep track of a visited array to skip the reverse edge.
    pub fn euler_path(&self, u: usize) -> Vec<usize> {
        let mut adj_iters = (0..self.num_v()).map(|u| self.adj_list(u))
                            .collect::<Vec<_>>();
        let mut edges = Vec::with_capacity(self.num_e());
        self.euler_recurse(u, &mut adj_iters, &mut edges);
        edges.reverse();
        edges
    }
}

// An iterator for convenient adjacency list traversal.
pub struct AdjListIterator<'a> {
    graph: &'a Graph,
    next_e: Option<usize>
}

impl<'a> Iterator for AdjListIterator<'a> {
    type Item = (usize, usize);
    
    // Produces an outgoing edge and vertex.
    fn next(&mut self) -> Option<Self::Item> {
        self.next_e.map( |e| {
            let v = self.graph.endp[e];
            self.next_e = self.graph.next[e];
            (e, v)
        })
    }
}

// Kruskal's minimum spanning tree algorithm on an undirected graph.
pub fn min_spanning_tree(graph: &Graph, weights: &[i64]) -> Vec<usize> {
    assert_eq!(graph.num_e(), 2 * weights.len());
    let mut edges = (0..weights.len()).collect::<Vec<_>>();
    edges.sort_by_key(|&i| weights[i]);
    
    let mut components = DisjointSets::new(graph.num_v());
    edges.into_iter()
        .filter(|&e| components.merge(graph.endp[2*e], graph.endp[2*e+1]))
        .collect()
}

#[derive(Clone)]
struct FlowVertex {
    lev: Option<usize>,
    cur: Option<usize> // TODO: consider making this an AdjListIterator
}

pub struct FlowEdge {
    pub flow: i64,
    pub cap: i64,
    pub cost: i64
}

pub struct FlowGraph {
    pub graph: Graph,
    vdata: Vec<FlowVertex>,
    pub edata: Vec<FlowEdge>,
}

impl FlowGraph {
    pub fn new(vmax: usize, emax: usize) -> FlowGraph {
        let data = FlowVertex { lev: None, cur: None };
        FlowGraph {
            graph: Graph::new(vmax, 2 * emax),
            vdata: vec![data; vmax],
            edata: Vec::with_capacity(2 * emax)
        }
    }
    
    pub fn add_edge(&mut self, a: usize, b: usize, cap: i64, cost: i64) {
        let data = FlowEdge { flow: 0, cap: cap, cost: cost };
        let rdata = FlowEdge { cost: -cost, ..data };
        self.edata.push(data);
        self.edata.push(rdata);
        self.graph.add_undirected_edge(a, b);
    }
    
    fn dfs(&mut self, u: usize, t: usize, f: i64) -> i64 {
        if u == t { return f; }
        let mut df = 0;
        
        while let Some(e) = self.vdata[u].cur {
            let v = self.graph.endp[e];
            if let (Some(lu), Some(lv)) = (self.vdata[u].lev, self.vdata[v].lev) {
                let rem_cap = self.edata[e].cap - self.edata[e].flow;
                if rem_cap > 0 && lv == lu + 1 {
                    let cf = self.dfs(v, t, ::std::cmp::min(rem_cap, f - df));
                    self.edata[e].flow += cf;
                    self.edata[e ^ 1].flow -= cf;
                    df += cf;
                    if df == f { break; }
                }
            }
            self.vdata[u].cur = self.graph.next[e];
        }
        return df;
    }
    
    fn bfs(&mut self, s: usize, t: usize) -> bool {
        for v in &mut self.vdata { v.lev = None; }
        let mut q = ::std::collections::VecDeque::<usize>::new();
        q.push_back(s);
        self.vdata[s].lev = Some(0);
        while let Some(u) = q.pop_front() {
            self.vdata[u].cur = self.graph.first[u];
            for (e, v) in self.graph.adj_list(u) {
                if self.vdata[v].lev == None && self.edata[e].flow < self.edata[e].cap {
                    q.push_back(v);
                    self.vdata[v].lev = Some(self.vdata[u].lev.unwrap() + 1);
                }
            }
        }
        self.vdata[t].lev != None
    }
    
    // Dinic's maximum flow / Hopcroft-Karp maximum bipartite matching: V^2E in
    // general, min(V^(2/3),sqrt(E))E on unit capacity, sqrt(V)E on bipartite.
    pub fn dinic(&mut self, s: usize, t: usize) -> i64 {
        let mut flow = 0;
        while self.bfs(s, t) {
            flow += self.dfs(s, t, 0x3f3f3f3f);
        }
        flow
    }
    
    pub fn min_cut(&self) -> Vec<usize> {
        (0..self.graph.num_e()).filter( |&e| {
            let u = self.graph.endp[e ^ 1];
            let v = self.graph.endp[e];
            self.vdata[u].lev.is_some() && self.vdata[v].lev.is_none()
        }).collect()
    }
}

// Strongly connected, 2-vertex-connected, and 2-edge-connected components
// should handle multiple-edges and self-loops
// USAGE: 1) new(); 2) add_edge(...); 3) compute_bcc();
// 4) use is_cut_vertex(vertex_index) or is_cut_edge(2 * edge_index)

#[derive(Clone)]
pub struct CCVertex {
    pub cc: usize,
    low: usize,
    vis: usize,
}

pub struct CCGraph<'a> {
    pub graph: &'a Graph,
    pub vdata: Vec<CCVertex>,
    pub vcc: Vec<usize>,
    pub n_cc: usize,
    pub n_vcc: usize,
    t: usize,
    verts: ::std::collections::VecDeque<usize>,
    //edges: ::std::collections::VecDeque<usize>
}

impl<'a> CCGraph<'a> {
    pub fn new(graph: &'a Graph) -> CCGraph {
        let data = CCVertex { cc: 0, low: 0, vis: 0 };
        let mut cc_graph = CCGraph {
            graph: graph,
            vdata: vec![data; graph.num_v()],
            vcc: vec![0; graph.num_e()],
            n_cc: 0,
            n_vcc: 0,
            t: 0,
            verts: ::std::collections::VecDeque::new(),
            //edges: ::std::collections::VecDeque::new()
        };
        for i in 0..graph.num_v() {
            if cc_graph.vdata[i].vis == 0 {
                cc_graph.scc(i);
            }
        }
        cc_graph
    }
    
    // SCCs form a DAG whose components are numbered in reverse topological order.
    fn scc(&mut self, u: usize) {
        self.t += 1;
        self.vdata[u].low = self.t;
        self.vdata[u].vis = self.t;
        self.verts.push_back(u);
        for (_, v) in self.graph.adj_list(u) {
            if self.vdata[v].vis == 0 { self.scc(v); }
            if self.vdata[v].cc == 0 {
                self.vdata[u].low = ::std::cmp::min(self.vdata[u].low, self.vdata[v].low);
            }
        }
        if self.vdata[u].vis <= self.vdata[u].low {
            self.n_cc += 1;
            while let Some(v) = self.verts.pop_back() {
                self.vdata[v].cc = self.n_cc;
                if v == u { break; }
            }
        }
    }
    
    pub fn two_sat_assign(&self) -> Option<Vec<bool>> {
        (0..self.graph.num_v()/2).map( |i| {
            let scc_true = self.vdata[2*i].cc;
            let scc_false = self.vdata[2*i+1].cc;
            if scc_true == scc_false { None } else { Some(scc_true < scc_false) }
        }).collect()
    }
    
    /*fn bcc(&mut self, u: usize, par: usize) {
        self.t += 1;
        self.vdata[u].low = self.t;
        self.vdata[u].vis = self.t;
        self.verts.push_back(u);
        for (e, v) in self.graph.adj_list(u) {
          if self.vdata[v] == None {
              self.edges.push_back(e);
              self.bcc(v, e);
              self.vdata[u].low = ::std::cmp::min(self.vdata[u].low, self.vdata[v].low);
              if self.vdata[u].vis <= self.vdata[v].low { // u is a cut vertex unless it's a one-child root
                  do {
                      let E = self.edges.top();
                      self.edges.pop_back();
                      vcc[E] = self.n_vcc;
                      vcc[E^1] = self.n_vcc;
                  } while e != E && e != (E^1);
                  self.n_vcc += 1;
              }
          }
          else if self.vdata[v].vis < self.vdata[u].vis && e != (par^1) {
              self.vdata[u].low = ::std::cmp::min(self.vdata[u].low, self.vdata[v].vis);
              self.edges.push_back(e);
          }
          else if v == u { // e is a self-loop
              self.vcc[e] = self.n_vcc;
              self.vcc[e ^ 1] = self.n_vcc;
              self.n_vcc += 1;
          }
        } 
        if self.vdata[u].vis <= self.vdata[u].low { // par is a cut edge unless par==-1
            do { v = self.verts.top(); self.verts.pop_back(); self.vdata[v].cc = self.n_cc; }
            while (v != u);
            self.n_cc += 1; 
        }
    }
    
    pub fn is_cut_vertex(&self, u: usize) -> bool {
        let vcc = self.vcc[self.graph.first[u]];
        for (e, _) in self.graph.adj_list(u) {
            if (self.vcc[e] != vcc) { return true; }
        }
        false
    }
    
    pub fn is_cut_edge(&self, e: usize) -> bool {
        let u = self.graph.endp[e ^ 1];
        let v = self.graph.endp[e];
        self.vdata[u].cc != self.vdata[v].cc
    }*/
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test_euler()
    {
        let mut graph = Graph::new(3, 4);
        graph.add_edge(0, 1);
        graph.add_edge(1, 0);
        graph.add_edge(1, 2);
        graph.add_edge(2, 1);
        assert_eq!(graph.euler_path(0), vec![0, 2, 3, 1]);
    }
    
    #[test]
    fn test_min_spanning_tree()
    {
        let mut graph = Graph::new(3, 3);
        graph.add_undirected_edge(0, 1);
        graph.add_undirected_edge(1, 2);
        graph.add_undirected_edge(2, 0);
        let weights = [7, 3, 5];
        let mst = min_spanning_tree(&graph, &weights);
        assert_eq!(mst, vec![1, 2]);
    }
    
    #[test]
    fn test_basic_flow()
    {
        let mut graph = FlowGraph::new(3, 2);
        graph.add_edge(0, 1, 4, 1);
        graph.add_edge(1, 2, 3, 1);
        let flow = graph.dinic(0, 2);
        assert_eq!(flow, 3);
    }
    
    #[test]
    fn test_two_sat()
    {
        let mut graph = Graph::new(6, 8);
        let (x, y, z) = (0, 2, 4);
        
        graph.add_two_sat_clause(x, z);
        graph.add_two_sat_clause(y^1, z^1);
        graph.add_two_sat_clause(y, y);
        assert_eq!(CCGraph::new(&graph).two_sat_assign(),
                   Some(vec![true, true, false]));
            
        graph.add_two_sat_clause(z, z);
        assert_eq!(CCGraph::new(&graph).two_sat_assign(), None);
    }
}
