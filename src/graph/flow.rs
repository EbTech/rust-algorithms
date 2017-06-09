use super::*;
use ::std::cmp::min;
const INF: i64 = 0x3f3f3f3f;

#[derive(Clone)]
struct FlowVertex {
    pot: i64,
    cur: Option<usize> // TODO: consider making this an AdjListIterator
}

pub struct FlowEdge {
    pub flow: i64, // TODO: should be the only mutable member while augmenting
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
        let data = FlowVertex { pot: 0, cur: None };
        FlowGraph {
            graph: Graph::new(vmax, 2 * emax),
            vdata: vec![data; vmax],
            edata: Vec::with_capacity(2 * emax)
        }
    }
    
    pub fn add_edge(&mut self, a: usize, b: usize, cap: i64, cost: i64) {
        let fdata = FlowEdge { flow: 0, cap: cap, cost: cost };
        let rdata = FlowEdge { flow: 0, cap: 0, cost: -cost };
        self.edata.push(fdata);
        self.edata.push(rdata);
        self.graph.add_undirected_edge(a, b);
    }
    
    // Dinic's maximum flow / Hopcroft-Karp maximum bipartite matching: V^2E in
    // general, min(V^(2/3),sqrt(E))E on unit capacity, sqrt(V)E on bipartite.
    pub fn dinic(&mut self, s: usize, t: usize) -> i64 {
        let mut flow = 0;
        while self.bfs(s, t) {
            flow += self.dfs(s, t, INF);
        }
        flow
    }
    
    fn bfs(&mut self, s: usize, t: usize) -> bool {
        for v in &mut self.vdata { v.pot = INF; }
        let mut q = ::std::collections::VecDeque::<usize>::new();
        q.push_back(s);
        self.vdata[s].pot = 0;
        while let Some(u) = q.pop_front() {
            self.vdata[u].cur = self.graph.first[u];
            for (e, v) in self.graph.adj_list(u) {
                if self.vdata[v].pot == INF && self.edata[e].flow < self.edata[e].cap {
                    q.push_back(v);
                    self.vdata[v].pot = self.vdata[u].pot + 1;
                }
            }
        }
        self.vdata[t].pot < INF
    }
    
    fn dfs(&mut self, u: usize, t: usize, f: i64) -> i64 {
        if u == t { return f; }
        let mut df = 0;
        
        while let Some(e) = self.vdata[u].cur {
            let v = self.graph.endp[e];
            let rem_cap = self.edata[e].cap - self.edata[e].flow;
            if rem_cap > 0 && self.vdata[v].pot == self.vdata[u].pot + 1 {
                let cf = self.dfs(v, t, min(rem_cap, f - df));
                self.edata[e].flow += cf;
                self.edata[e ^ 1].flow -= cf;
                df += cf;
                if df == f { break; }
            }
            self.vdata[u].cur = self.graph.next[e];
        }
        return df;
    }
    
    // After running maximum flow, use this to recover the dual minimum cut.
    pub fn min_cut(&self) -> Vec<usize> {
        (0..self.graph.num_e()).filter( |&e| {
            let u = self.graph.endp[e ^ 1];
            let v = self.graph.endp[e];
            self.vdata[u].pot < INF && self.vdata[v].pot == INF
        }).collect()
    }
    
    // Minimum cost maximum flow, assuming no negative-cost cycles.
    pub fn mcf(&mut self, s: usize, t: usize) -> (i64, i64) {
        for u in 0..self.graph.num_v() { self.vdata[u].pot = 0; }
        
        // If there are no negative-cost edges, this Bellman-Ford can be omitted.
        for _ in 1..self.graph.num_v() {
            for e in 1..self.graph.num_e() {
                if self.edata[e].cap > 0 {
                    let u = self.graph.endp[e ^ 1];
                    let v = self.graph.endp[e];
                    self.vdata[v].pot = min(self.vdata[v].pot, self.vdata[u].pot + self.edata[e].cost);
                }
            }
        }
        
        let (mut cost, mut flow) = (0, 0);
        while let Some((dc, df)) = self.mcf_augment(s, t) {
            cost += dc;
            flow += df;
        }
        (cost, flow)
    }
    
    fn mcf_augment(&mut self, s: usize, t: usize) -> Option<(i64, i64)> {
        let mut vis = vec![false; self.graph.num_v()];
        let mut dist = vec![INF; self.graph.num_v()];
        let mut par = vec![None; self.graph.num_v()];
        dist[s] = 0; // Do Dijkstra.
        while let Some(u) = (0..self.graph.num_v()).filter(|&u| !vis[u])
                            .min_by_key(|&u| dist[u]) {
            vis[u] = true;
            for (e, v) in self.graph.adj_list(u) {
                if self.edata[e].flow < self.edata[e].cap {
                    let d = dist[u] + self.vdata[u].pot - self.vdata[v].pot + self.edata[e].cost;
                    if dist[v] > d {
                        dist[v] = d;
                        par[v] = Some(e);
                    }
                }
            }
        }
        if dist[t] >= INF { return None; }
        let (mut dc, mut df) = (0, INF);
        let mut i = t;
        while let Some(e) = par[i] {
            df = min(df, self.edata[e].cap - self.edata[e].flow);
            i = self.graph.endp[e ^ 1];
        }
        i = t;
        while let Some(e) = par[i] {
            self.edata[e].flow += df;
            self.edata[e ^ 1].flow -= df;
            dc += df * self.edata[e].cost;
            i = self.graph.endp[e ^ 1];
        }
        for u in 0..self.graph.num_v() {
            self.vdata[u].pot = min(INF, dist[u] + self.vdata[u].pot);
        }
        Some((df, dc))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    
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
    fn test_min_cost_flow()
    {
        let mut graph = FlowGraph::new(4, 4);
        graph.add_edge(0, 1, 10, -10);
        graph.add_edge(1, 2, 7, 8);
        graph.add_edge(2, 3, 7, 8);
        graph.add_edge(1, 3, 7, 10);
        let (flow, cost) = graph.mcf(0, 3);
        assert_eq!(flow, 10);
        assert_eq!(cost, 18);
    }
}
