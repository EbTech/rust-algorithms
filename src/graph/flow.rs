//! Maximum flows, matchings, and minimum cuts.
use super::{AdjListIterator, Graph};

/// Representation of a network flow problem with (optional) costs.
pub struct FlowGraph {
    /// Owned graph, managed by this FlowGraph object.
    pub graph: Graph,
    /// Edge capacities.
    pub cap: Vec<i64>,
    /// Edge cost per unit flow.
    pub cost: Vec<i64>,
}

impl FlowGraph {
    /// An upper limit to the flow.
    const INF: i64 = 0x3f3f_3f3f_3f3f_3f3f;

    /// Initializes an flow network with vmax vertices and no edges.
    pub fn new(vmax: usize, emax_hint: usize) -> Self {
        Self {
            graph: Graph::new(vmax, 2 * emax_hint),
            cap: Vec::with_capacity(2 * emax_hint),
            cost: Vec::with_capacity(2 * emax_hint),
        }
    }

    /// Adds an edge with specified capacity and cost. The reverse edge is also
    /// added for residual graph computation, but has zero capacity.
    pub fn add_edge(&mut self, u: usize, v: usize, cap: i64, cost: i64) {
        self.cap.push(cap);
        self.cap.push(0);
        self.cost.push(cost);
        self.cost.push(-cost);
        self.graph.add_undirected_edge(u, v);
    }

    /// Dinic's algorithm to find the maximum flow from s to t where s != t.
    /// Generalizes the Hopcroft-Karp maximum bipartite matching algorithm.
    /// V^2E in general, min(V^(2/3),sqrt(E))E when all edges are unit capacity,
    /// sqrt(V)E when all vertices are unit capacity as in bipartite graphs.
    ///
    /// # Panics
    ///
    /// Panics if the maximum flow is 2^63 or larger.
    pub fn dinic(&self, s: usize, t: usize) -> (i64, Vec<i64>) {
        let mut flow = vec![0; self.graph.num_e()];
        let mut max_flow = 0;
        loop {
            let dist = self.dinic_search(s, &flow);
            if dist[t] == Self::INF {
                break;
            }
            // Keep track of adjacency lists to avoid revisiting blocked edges.
            let mut adj_iters = (0..self.graph.num_v())
                .map(|u| self.graph.adj_list(u).peekable())
                .collect::<Vec<_>>();
            max_flow += self.dinic_augment(s, t, Self::INF, &dist, &mut adj_iters, &mut flow);
        }
        (max_flow, flow)
    }

    // Compute BFS distances to restrict attention to shortest path edges.
    fn dinic_search(&self, s: usize, flow: &[i64]) -> Vec<i64> {
        let mut dist = vec![Self::INF; self.graph.num_v()];
        let mut q = ::std::collections::VecDeque::new();
        dist[s] = 0;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            for (e, v) in self.graph.adj_list(u) {
                if dist[v] == Self::INF && flow[e] < self.cap[e] {
                    dist[v] = dist[u] + 1;
                    q.push_back(v);
                }
            }
        }
        dist
    }

    // Pushes a blocking flow that increases the residual's s-t distance.
    fn dinic_augment(
        &self,
        u: usize,
        t: usize,
        f: i64,
        dist: &[i64],
        adj: &mut [::std::iter::Peekable<AdjListIterator>],
        flow: &mut [i64],
    ) -> i64 {
        if u == t {
            return f;
        }
        let mut df = 0;

        while let Some(&(e, v)) = adj[u].peek() {
            let rem_cap = (self.cap[e] - flow[e]).min(f - df);
            if rem_cap > 0 && dist[v] == dist[u] + 1 {
                let cf = self.dinic_augment(v, t, rem_cap, dist, adj, flow);
                flow[e] += cf;
                flow[e ^ 1] -= cf;
                df += cf;
                if df == f {
                    break;
                }
            }
            // The current edge is either saturated or blocked.
            adj[u].next();
        }
        df
    }

    /// After running maximum flow, use this to recover the dual minimum cut.
    pub fn min_cut(&self, dist: &[i64]) -> Vec<usize> {
        (0..self.graph.num_e())
            .filter(|&e| {
                let u = self.graph.endp[e ^ 1];
                let v = self.graph.endp[e];
                dist[u] < Self::INF && dist[v] == Self::INF
            })
            .collect()
    }

    /// Among all s-t maximum flows, finds one with minimum cost, assuming
    /// s != t and no negative-cost cycles.
    ///
    /// # Panics
    ///
    /// Panics if the flow or cost overflow a 64-bit signed integer.
    pub fn mcf(&self, s: usize, t: usize) -> (i64, i64, Vec<i64>) {
        let mut pot = vec![0; self.graph.num_v()];

        // Bellman-Ford deals with negative-cost edges at initialization.
        for _ in 1..self.graph.num_v() {
            for e in 0..self.graph.num_e() {
                if self.cap[e] > 0 {
                    let u = self.graph.endp[e ^ 1];
                    let v = self.graph.endp[e];
                    pot[v] = pot[v].min(pot[u] + self.cost[e]);
                }
            }
        }

        let mut flow = vec![0; self.graph.num_e()];
        let (mut min_cost, mut max_flow) = (0, 0);
        loop {
            let par = self.mcf_search(s, &flow, &mut pot);
            if par[t] == None {
                break;
            }
            let (dc, df) = self.mcf_augment(t, &par, &mut flow);
            min_cost += dc;
            max_flow += df;
        }
        (min_cost, max_flow, flow)
    }

    // Maintains Johnson's potentials to prevent negative-cost residual edges.
    // This allows running Dijkstra instead of the slower Bellman-Ford.
    fn mcf_search(&self, s: usize, flow: &[i64], pot: &mut [i64]) -> Vec<Option<usize>> {
        let mut vis = vec![false; self.graph.num_v()];
        let mut dist = vec![Self::INF; self.graph.num_v()];
        let mut par = vec![None; self.graph.num_v()];

        dist[s] = 0;
        while let Some(u) = (0..self.graph.num_v())
            .filter(|&u| !vis[u] && dist[u] < Self::INF)
            .min_by_key(|&u| dist[u] - pot[u])
        {
            vis[u] = true;
            pot[u] = dist[u];
            for (e, v) in self.graph.adj_list(u) {
                if dist[v] > dist[u] + self.cost[e] && flow[e] < self.cap[e] {
                    dist[v] = dist[u] + self.cost[e];
                    par[v] = Some(e);
                }
            }
        }
        par
    }

    // Pushes flow along an augmenting path of minimum cost.
    fn mcf_augment(&self, t: usize, par: &[Option<usize>], flow: &mut [i64]) -> (i64, i64) {
        let (mut dc, mut df) = (0, Self::INF);
        let mut u = t;
        while let Some(e) = par[u] {
            df = df.min(self.cap[e] - flow[e]);
            u = self.graph.endp[e ^ 1];
        }
        u = t;
        while let Some(e) = par[u] {
            flow[e] += df;
            flow[e ^ 1] -= df;
            dc += df * self.cost[e];
            u = self.graph.endp[e ^ 1];
        }
        (dc, df)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_basic_flow() {
        let mut graph = FlowGraph::new(3, 2);
        graph.add_edge(0, 1, 4, 1);
        graph.add_edge(1, 2, 3, 1);

        let flow = graph.dinic(0, 2).0;
        assert_eq!(flow, 3);
    }

    #[test]
    fn test_min_cost_flow() {
        let mut graph = FlowGraph::new(4, 4);
        graph.add_edge(0, 1, 10, -10);
        graph.add_edge(1, 2, 7, 8);
        graph.add_edge(2, 3, 7, 8);
        graph.add_edge(1, 3, 7, 10);

        let (cost, flow, _) = graph.mcf(0, 3);
        assert_eq!(cost, 18);
        assert_eq!(flow, 10);
    }

    #[test]
    fn test_max_matching() {
        let mut graph = FlowGraph::new(14, 4);

        let source = 0;
        let sink = 13;

        //Vertex indices of "left hand side" of bipartite graph go from [left_start, right_start)
        let left_start = 1;
        //Vertex indices of "right hand side" of bipartite graph go from [right_start, sink)
        let right_start = 7;

        //Initialize source / sink connections; both left & right have 6 nodes
        for lhs_vertex in left_start..left_start + 6 {
            graph.add_edge(source, lhs_vertex, 1, 1);
        }

        for rhs_vertex in right_start..right_start + 6 {
            graph.add_edge(rhs_vertex, sink, 1, 1);
        }

        graph.add_edge(left_start + 0, right_start + 1, 1, 1);
        graph.add_edge(left_start + 0, right_start + 2, 1, 1);
        graph.add_edge(left_start + 2, right_start + 0, 1, 1);
        graph.add_edge(left_start + 2, right_start + 3, 1, 1);
        graph.add_edge(left_start + 3, right_start + 2, 1, 1);
        graph.add_edge(left_start + 4, right_start + 2, 1, 1);
        graph.add_edge(left_start + 4, right_start + 3, 1, 1);
        graph.add_edge(left_start + 5, right_start + 5, 1, 1);

        let (flow_amt, flow) = graph.dinic(source, sink);
        assert_eq!(flow_amt, 5);

        //L->R edges in maximum matching
        let left_right_edges = flow
            .into_iter()
            .enumerate()
            .filter(|&(_e, f)| f > 0)
            //map to u->v
            .map(|(e, _f)| (graph.graph.endp[e ^ 1], graph.graph.endp[e]))
            //leave out source and sink nodes
            .filter(|&(u, v)| u != source && v != sink)
            .collect::<Vec<_>>();

        assert_eq!(
            left_right_edges,
            vec![(1, 8), (3, 7), (4, 9), (5, 10), (6, 12)]
        );
    }
}
