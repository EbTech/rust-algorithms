//use std::cmp::min;
//const INF: i32 = 0x3f3f3f3f;

pub struct DisjointSets {
    parent: Vec<usize>
}

impl DisjointSets {
    pub fn new(size: usize) -> DisjointSets {
        DisjointSets { parent: (0..size).collect() }
    }
    
    pub fn find(&mut self, u: usize) -> usize {
        let par_u = self.parent[u];
        if par_u != u { self.parent[u] = self.find(par_u); }
        self.parent[u]
    }
    
    pub fn union(&mut self, u: usize, v: usize) {
        let (pu, pv) = (self.find(u), self.find(v));
        if pu != pv { self.parent[pu] = pv; }
    }
}

#[derive(Clone)]
pub struct Vertex<V> {
    pub first: Option<usize>,
    pub data: V
}

pub struct Edge<E> {
    pub endp: usize,
    pub next: Option<usize>,
    pub data: E
}

pub struct Graph<V, E> {
    pub verts: Vec<Vertex<V>>,
    pub edges: Vec<Edge<E>>,
}

impl<V, E> Graph<V, E> where V: Clone {
    pub fn new(vmax: usize, emax: usize, vdata: V) -> Graph<V, E> {
        Graph {
            verts: vec![Vertex{ first: None, data: vdata }; vmax],
            edges: Vec::with_capacity(2 * emax)
        }
    }
    
    pub fn add_edge(&mut self, a: usize, b: usize, data: E, rdata: E) {
        self.edges.push(Edge{
            endp: a,
            next: self.verts[a].first,
            data: data
        });
        self.verts[a].first = Some(self.edges.len() - 1);
        self.edges.push(Edge{
            endp: b,
            next: self.verts[b].first,
            data: rdata
        });
        self.verts[b].first = Some(self.edges.len() - 1);
    }
}

#[derive(Clone)]
struct VFlowData {
    lev: Option<usize>,
    cur: Option<usize>
}

pub struct EFlowData {
    pub flow: i64,
    pub cap: i64,
    pub cost: i64
}

pub struct FlowGraph {
    graph: Graph<VFlowData, EFlowData>
}

impl FlowGraph {
    pub fn new(vmax: usize, emax: usize) -> FlowGraph {
        let data = VFlowData { lev: None, cur: None };
        let graph = Graph::new(vmax, emax, data);
        FlowGraph { graph: graph }
    }
    
    pub fn add_edge(&mut self, a: usize, b: usize, cap: i64, cost: i64) {
        let data = EFlowData { flow: 0, cap: cap, cost: cost };
        let rdata = EFlowData { cost: -cost, ..data };
        self.graph.add_edge(a, b, data, rdata);
    }
    
    fn bfs(&mut self, s: usize, t: usize) -> bool {
        let g = &mut self.graph;
        for v in &mut g.verts { v.data.lev = None; }
        let mut q = ::std::collections::VecDeque::<usize>::new();
        q.push_back(s);
        g.verts[s].data.lev = Some(0);
        while let Some(u) = q.pop_front() {
            g.verts[u].data.cur = g.verts[u].first;
            let mut e = g.verts[u].first;
            while let Some(edge_id) = e {
                let edge = &g.edges[edge_id];
                let rev_edge = &g.edges[edge_id ^ 1];
                let v = rev_edge.endp;
                if g.verts[v].data.lev == None && edge.data.flow < edge.data.cap {
                    q.push_back(v);
                    g.verts[v].data.lev = Some(g.verts[u].data.lev.unwrap() + 1);
                }
                e = g.edges[edge_id].next;
            }
        }
        return g.verts[t].data.lev != None;
    }
    
    fn dfs(&mut self, u: usize, t: usize, f: i64) -> i64 {
        if u == t { return f; }
        let mut df = 0;
        
        while let Some(edge_id) = self.graph.verts[u].data.cur {
            let v = self.graph.edges[edge_id ^ 1].endp;
            if let (Some(lu), Some(lv)) = (self.graph.verts[u].data.lev, self.graph.verts[v].data.lev) {
                let rem_cap = self.graph.edges[edge_id].data.cap
                            - self.graph.edges[edge_id].data.flow;
                if rem_cap > 0 && lv == lu + 1 {
                    let cf = self.dfs(v, t, ::std::cmp::min(rem_cap, f - df));
                    self.graph.edges[edge_id].data.flow += cf;
                    self.graph.edges[edge_id ^ 1].data.flow -= cf;
                    df += cf;
                    if df == f { break; }
                }
            }
            self.graph.verts[u].data.cur = self.graph.edges[edge_id].next;
        }
        return df;
    }
    
    // Dinic's fast maximum flow: V^2E in general,
    // min(V^(2/3),sqrt(E))E on unit caps, sqrt(V)E on bipartite
    pub fn dinic(&mut self, s: usize, t: usize) -> i64 {
        let mut flow = 0;
        while self.bfs(s, t) {
            flow += self.dfs(s, t, 0x3f3f3f3f);
        }
        flow
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
}
