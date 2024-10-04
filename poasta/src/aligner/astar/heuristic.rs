use super::AlignableGraph;


pub trait AstarHeuristic<G, T>: Default
where
    G: AlignableGraph,
{
    fn init(&mut self, graph: &G, seq: &[u8]);
    
    fn h(&self, item: &T) -> usize;
}

#[derive(Default)]
pub struct Dijkstra;

impl<G, T> AstarHeuristic<G, T> for Dijkstra
where 
    G: AlignableGraph
{
    #[inline(always)]
    fn init(&mut self, _: &G, _: &[u8]) {
        
    }
    
    #[inline(always)]
    fn h(&self, _: &T) -> usize {
        0
    }
}


#[derive(Default, Clone)]
pub struct MinGapCostAffine {
    min_dist_to_end: Vec<usize>,
    max_dist_to_end: Vec<usize>,
}

impl<G, T> AstarHeuristic<G, T> for MinGapCostAffine
where
    G: AlignableGraph
{
    fn init(&mut self, graph: &G, _: &[u8]) {
        self.min_dist_to_end = vec![0; graph.node_count()];
        self.max_dist_to_end = vec![0; graph.node_count()];
    }
    
    fn h(&self, _: &T) -> usize {
        0
    }
}