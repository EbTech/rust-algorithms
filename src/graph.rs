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

pub struct Graph {
    pub first: Vec<Option<usize>>,
    pub next: Vec<Option<usize>>,
    pub endp: Vec<usize>,
}

impl Graph {
    pub fn new(vmax: usize, emax: usize) -> Graph {
        Graph {
            first: vec![None; vmax],
            next: Vec::with_capacity(2 * emax),
            endp: Vec::with_capacity(2 * emax)
        }
    }
    
    pub fn add_edge(&mut self, a: usize, b: usize) {
        for &u in &[a, b] {
            self.next.push(self.first[u]);
            self.first[u] = Some(self.endp.len());
            self.endp.push(u);
        }
    }
    
    pub fn adj_list<'a>(&'a self, u: usize) -> AdjListIterator<'a> {
        AdjListIterator { graph: self, next_e: self.first[u] }
    }
}

pub struct AdjListIterator<'a> {
    graph: &'a Graph,
    next_e: Option<usize>
}

impl<'a> ::std::iter::Iterator for AdjListIterator<'a> {
    // produces an outgoing edge and vertex
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        self.next_e.map( |e| {
            let v = self.graph.endp[e ^ 1];
            self.next_e = self.graph.next[e];
            (e, v)
        })
    }
}

#[derive(Clone)]
pub struct FlowVertex {
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
    pub vdata: Vec<FlowVertex>,
    pub edata: Vec<FlowEdge>,
}

impl FlowGraph {
    pub fn new(vmax: usize, emax: usize) -> FlowGraph {
        let data = FlowVertex { lev: None, cur: None };
        FlowGraph {
            graph: Graph::new(vmax, emax),
            vdata: vec![data; vmax],
            edata: Vec::with_capacity(2 * emax)
        }
    }
    
    pub fn add_edge(&mut self, a: usize, b: usize, cap: i64, cost: i64) {
        let data = FlowEdge { flow: 0, cap: cap, cost: cost };
        let rdata = FlowEdge { cost: -cost, ..data };
        self.edata.push(data);
        self.edata.push(rdata);
        self.graph.add_edge(a, b);
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
        return self.vdata[t].lev != None;
    }
    
    fn dfs(&mut self, u: usize, t: usize, f: i64) -> i64 {
        if u == t { return f; }
        let mut df = 0;
        
        while let Some(e) = self.vdata[u].cur {
            let v = self.graph.endp[e ^ 1];
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
