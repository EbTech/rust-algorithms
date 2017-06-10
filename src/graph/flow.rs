use graph::Graph;
use graph::AdjListIterator;
use ::std::cmp::min;
const INF: i64 = 0x3f3f3f3f;

// Representation of a network flow problem with (optional) costs.
pub struct FlowGraph {
    pub graph: Graph,
    pub cap: Vec<i64>,
    pub cost: Vec<i64>
}

impl FlowGraph {
    // Initializes an flow network with vmax vertices and no edges.
    pub fn new(vmax: usize, emax: usize) -> FlowGraph {
        FlowGraph {
            graph: Graph::new(vmax, 2 * emax),
            cap: Vec::with_capacity(2 * emax),
            cost: Vec::with_capacity(2 * emax)
        }
    }
    
    // Adds an edge with specified capacity and cost. The reverse edge is also
    // added for residual graph computation, but has zero capacity.
    pub fn add_edge(&mut self, a: usize, b: usize, cap: i64, cost: i64) {
        self.cap.push(cap); self.cost.push(cost);
        self.cap.push(0); self.cost.push(-cost);
        self.graph.add_undirected_edge(a, b);
    }
    
    // Dinic's maximum flow / Hopcroft-Karp maximum bipartite matching: V^2E in
    // general, min(V^(2/3),sqrt(E))E on unit capacity, sqrt(V)E on bipartite.
    pub fn dinic(&self, s: usize, t: usize) -> i64 {
        let mut flow = vec![0; self.graph.num_e()];
        let mut max_flow = 0;
        while let Some(df) = self.dinic_augment(s, t, &mut flow) {
            max_flow += df;
        }
        max_flow
    }
    
    // Pushes a saturating flow that increases the residual's s-t distance.
    pub fn dinic_augment(&self, s: usize, t: usize, flow: &mut [i64]) -> Option<i64> {
        let mut dist = vec![INF; self.graph.num_v()];
        let mut q = ::std::collections::VecDeque::new();
        dist[s] = 0;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            for (e, v) in self.graph.adj_list(u) {
                if dist[v] == INF && flow[e] < self.cap[e] {
                    dist[v] = dist[u] + 1;
                    q.push_back(v);
                }
            }
        }
        if dist[t] == INF { return None; }
        
        let mut adj_iters = (0..self.graph.num_v())
                            .map(|u| self.graph.adj_list(u))
                            .collect::<Vec<_>>();
        let df = self.dinic_dfs(s, t, INF, &dist, &mut adj_iters, flow);
        Some(df)
    }
    
    fn dinic_dfs(&self, u: usize, t: usize, f: i64, dist: &[i64],
                 adj: &mut [AdjListIterator], flow: &mut [i64]) -> i64 {
        if u == t { return f; }
        let mut df = 0;
        
        while let Some((e, v)) = adj[u].next() {
            let rem_cap = min(self.cap[e] - flow[e], f - df);
            if rem_cap > 0 && dist[v] == dist[u] + 1 {
                let cf = self.dinic_dfs(v, t, rem_cap, dist, adj, flow);
                flow[e] += cf;
                flow[e ^ 1] -= cf;
                df += cf;
                if df == f { break; }
            }
        }
        return df;
    }
    
    // After running maximum flow, use this to recover the dual minimum cut.
    pub fn min_cut(&self, dist: &[i64]) -> Vec<usize> {
        (0..self.graph.num_e()).filter( |&e| {
            let u = self.graph.endp[e ^ 1];
            let v = self.graph.endp[e];
            dist[u] < INF && dist[v] == INF
        }).collect()
    }
    
    // Minimum cost maximum flow, assuming no negative-cost cycles.
    pub fn mcf(&self, s: usize, t: usize) -> (i64, i64) {
        let mut pot = vec![0; self.graph.num_v()];
        
        // Bellman-Ford needed here to deal with negative-cost edges.
        for _ in 1..self.graph.num_v() {
            for e in 0..self.graph.num_e() {
                if self.cap[e] > 0 {
                    let u = self.graph.endp[e ^ 1];
                    let v = self.graph.endp[e];
                    pot[v] = min(pot[v], pot[u] + self.cost[e]);
                }
            }
        }
        
        let mut flow = vec![0; self.graph.num_e()];
        let (mut min_cost, mut max_flow) = (0, 0);
        while let Some((dc, df)) = self.mcf_augment(s, t, &mut pot, &mut flow) {
            min_cost += dc;
            max_flow += df;
        }
        (min_cost, max_flow)
    }
    
    // Pushes along an augmenting path of minimum cost while maintaining the
    // vertex potentials so that no negative-weight residual edges appear.
    pub fn mcf_augment(&self, s: usize, t: usize, pot: &mut [i64], flow: &mut [i64]) -> Option<(i64, i64)> {
        let mut vis = vec![false; self.graph.num_v()];
        let mut dist = vec![INF; self.graph.num_v()];
        let mut par = vec![None; self.graph.num_v()];
        
        // Potential-reweighted graph has no negative edges, so run Dijkstra.
        dist[s] = 0;
        while let Some(u) = (0..self.graph.num_v()).filter(|&u| !vis[u])
                            .min_by_key(|&u| dist[u]) {
            vis[u] = true;
            for (e, v) in self.graph.adj_list(u) {
                let d = dist[u] + pot[u] - pot[v] + self.cost[e];
                if dist[v] > d && flow[e] < self.cap[e] {
                    dist[v] = d;
                    par[v] = Some(e);
                }
            }
        }
        if dist[t] == INF { return None; }
        
        let (mut dc, mut df) = (0, INF);
        let mut u = t;
        while let Some(e) = par[u] {
            df = min(df, self.cap[e] - flow[e]);
            u = self.graph.endp[e ^ 1];
        }
        u = t;
        while let Some(e) = par[u] {
            flow[e] += df;
            flow[e ^ 1] -= df;
            dc += df * self.cost[e];
            u = self.graph.endp[e ^ 1];
        }
        for u in 0..self.graph.num_v() {
            pot[u] = min(INF, pot[u] + dist[u]);
        }
        Some((dc, df))
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
        let (cost, flow) = graph.mcf(0, 3);
        assert_eq!(cost, 18);
        assert_eq!(flow, 10);
    }
}
