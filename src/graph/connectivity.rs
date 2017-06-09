use super::*;
use ::std::cmp::min;

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
    verts: Vec<usize>,
    //edges: Vec<usize>
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
            verts: Vec::new(),
            //edges: Vec::new()
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
        self.verts.push(u);
        for (_, v) in self.graph.adj_list(u) {
            if self.vdata[v].vis == 0 { self.scc(v); }
            if self.vdata[v].cc == 0 {
                self.vdata[u].low = min(self.vdata[u].low, self.vdata[v].low);
            }
        }
        if self.vdata[u].vis <= self.vdata[u].low {
            self.n_cc += 1;
            while let Some(v) = self.verts.pop() {
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
    
    // Biconnected components are a work in progress.
    /*fn bcc(&mut self, u: usize, par: usize) {
        self.t += 1;
        self.vdata[u].low = self.t;
        self.vdata[u].vis = self.t;
        self.verts.push(u);
        for (e, v) in self.graph.adj_list(u) {
          if self.vdata[v] == None {
              self.edges.push_back(e);
              self.bcc(v, e);
              self.vdata[u].low = min(self.vdata[u].low, self.vdata[v].low);
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
              self.vdata[u].low = min(self.vdata[u].low, self.vdata[v].vis);
              self.edges.push_back(e);
          }
          else if v == u { // e is a self-loop
              self.vcc[e] = self.n_vcc;
              self.vcc[e ^ 1] = self.n_vcc;
              self.n_vcc += 1;
          }
        } 
        if self.vdata[u].vis <= self.vdata[u].low { // par is a cut edge unless par==-1
            while let Some(v) = self.verts.pop() {
                self.vdata[v].cc = self.n_cc;
                if v == u { break; }
            }
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
