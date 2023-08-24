use std::collections::VecDeque;
use serde::{Deserialize, Serialize};
use crate::aligner::extend::PathExtender;
use crate::graphs::{AlignableGraph, NodeIndexType};
use crate::aligner::offsets::{Diag, DiagonalPoint, DiagRange, MatrixPoint, OffsetType};
use crate::aligner::visited::{merge_intervals, VisitedInterval, VisitedIntervalType,
                              VisitedIntervalData, VisitedIntervalTree};


/// Holds the visited offsets for each diagonal using a set of interval trees.
#[derive(Default, Clone, Debug, Serialize, Deserialize)]
pub(crate) struct Layer<O>
where
    O: OffsetType,
{
    k_min: Diag,
    k_max: Diag,
    diagonals: VecDeque<VisitedIntervalTree<O>>
}

impl<O> Layer<O>
where
    O: OffsetType,
{
    pub fn new(mut k_min: Diag, mut k_max: Diag, mut diagonals: VecDeque<VisitedIntervalTree<O>>) -> Self {

        // Remove empty interval trees at the outer diagonals
        while let Some(front) = diagonals.front() {
            if front.is_empty() {
                k_min += 1;
                diagonals.pop_front();
            } else {
                break;
            }
        }

        while let Some(back) = diagonals.back() {
            if back.is_empty() {
                k_max -= 1;
                diagonals.pop_back();
            } else {
                break;
            }
        }

        Self { k_min, k_max, diagonals }
    }

    pub fn initial<N>(start_nodes: &[N]) -> Self
    where
        N: NodeIndexType,
    {
        if start_nodes.is_empty() {
            panic!("No start nodes given!")
        }

        let k_min = start_nodes.iter()
            .map(|n| (*n, O::zero()).diag())
            .min()
            .unwrap();

        let k_max = start_nodes.iter()
            .map(|n| (*n, O::zero()).diag())
            .max()
            .unwrap();

        eprintln!("INITIAL k_min: {:?}, k_max: {:?}", k_min, k_max);

        let num = k_max.diff(&k_min) + 1;
        let mut diagonals: VecDeque<VisitedIntervalTree<O>> = VecDeque::with_capacity(num);
        for _ in 0..num {
            diagonals.push_back(VisitedIntervalTree::new());
        }

        for node in start_nodes {
            let diag = (*node, O::zero()).diag();
            let ix = diag.diff(&k_min);
            diagonals[ix].insert(VisitedIntervalType::new_size_one(O::zero(), VisitedIntervalData::initial()));
        }

        for diag in &mut diagonals {
            diag.sort_and_index();
        }

        Self {
            k_min,
            k_max,
            diagonals,
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.diagonals.is_empty()
    }

    pub(crate) fn kmin(&self) -> Option<Diag> {
        if self.is_empty() { None } else { Some(self.k_min) }
    }

    pub(crate) fn kmax(&self) -> Option<Diag> {
        if self.is_empty() { None } else { Some(self.k_max) }
    }

    pub(crate) fn get_diagonals(&self) -> &VecDeque<VisitedIntervalTree<O>> {
        &self.diagonals
    }

    pub(crate) fn get_diagonals_mut(&mut self) -> &mut VecDeque<VisitedIntervalTree<O>> {
        &mut self.diagonals
    }

    /// Convert a diagonal to an index in our [`diagonals'] vector
    pub(crate) fn diag_to_ix(&self, diag: Diag) -> Option<usize> {
        if self.is_empty() {
            return None;
        }

        if diag < self.k_min || diag > self.k_max {
            None
        } else {
            Some(diag.diff(&self.k_min))
        }
    }

    /// Check current interval end points, and check their successor nodes in the graph. Identify
    /// the lowest diagonal of any of those successors.
    pub(crate) fn get_successors_kmin<G: AlignableGraph<NodeIndex=N>, N: NodeIndexType>(&self, graph: &G) -> Option<Diag> {
        if self.is_empty() {
            return None;
        }

        DiagRange::closed(self.k_min, self.k_max).into_iter().zip(&self.diagonals)
            .filter_map(|(diag, itree)| {
                itree.iter()
                    .filter_map(|ival| {
                        eprintln!("{:?}", ival);
                        let offset = ival.end - O::one();
                        let row: G::NodeIndex = (diag, offset).row();
                        let node = graph.row_to_node_ix(row);

                        graph.successors(node)
                            .map(|succ| (succ, ival.end).diag())
                            .min()
                    })
                    .min()
            })
            .min()
    }

    /// Check current interval end points, and check their successor nodes in the graph. Identify
    /// the highest diagonal of any of those successors.
    pub(crate) fn get_successors_kmax<G: AlignableGraph<NodeIndex=N>, N: NodeIndexType>(&self, graph: &G) -> Option<Diag> {
        if self.is_empty() {
            return None;
        }

        DiagRange::closed(self.k_min, self.k_max).into_iter().zip(&self.diagonals)
            .filter_map(|(diag, itree)| {
                itree.iter()
                    .filter_map(|ival| {
                        let offset = ival.end - O::one();
                        let row: G::NodeIndex = (diag, offset).row();
                        let node = graph.row_to_node_ix(row);

                        graph.successors(node)
                            .map(|succ| (succ, ival.end).diag())
                            .max()
                    })
                    .max()
            })
            .max()
    }

    pub(crate) fn get_visited(&self, diag: Diag) -> Option<&VisitedIntervalTree<O>> {
        if self.is_empty() {
            return None;
        }

        self.diag_to_ix(diag).map(|ix| &self.diagonals[ix])
    }

    pub(crate) fn get_visited_mut(&mut self, diag: Diag) -> Option<&mut VisitedIntervalTree<O>> {
        if self.is_empty() {
            return None;
        }

        self.diag_to_ix(diag).map(|ix| &mut self.diagonals[ix])
    }

    pub fn visited(&self, diag: Diag, offset: O) -> bool {
        if self.is_empty() {
            return false;
        }

        self.get_visited(diag)
            .map_or(false, |tree| tree.contains(offset))
    }

    pub(crate) fn extend<G>(&mut self, graph: &G, seq: &[u8]) -> Option<(Diag, &VisitedIntervalType<O>)>
    where
        G: AlignableGraph
    {
        if self.is_empty() {
            return None;
        }

        eprintln!("CURRENT: {:?}", self.diagonals);

        let mut extender = PathExtender::new(graph, seq, self);
        for (diag, itree) in DiagRange::closed(self.k_min, self.k_max)
            .into_iter()
            .zip(&self.diagonals)
        {
            eprintln!("Extending diagonal {:?}", diag);
            for ival in itree.iter() {
                if ival.extendable() {
                    extender.extend(diag, ival);
                }
            }
        }

        eprintln!("New extended intervals: {:?}", extender.get_new_intervals());

        // Merge new intervals with current
        if !extender.get_new_intervals().is_empty() {
            let mut new_ivals = extender.into_new_intervals();
            eprintln!("New intervals: {:?}", new_ivals);
            let new_kmin = *new_ivals.keys().min().unwrap();
            let new_kmax = *new_ivals.keys().max().unwrap();

            // First add intervals on new diagonals, not present in the layer before
            if new_kmin < self.k_min {
                let kmin_diff = self.k_min.diff(&new_kmin);
                eprintln!("k_min old: {:?}, k_min new: {:?}, diff: {:?}", self.k_min, new_kmin, kmin_diff);

                for i in 0..kmin_diff {
                    let new_diag = self.k_min - i as i64 - 1;
                    eprintln!("adding {:?}", new_diag);
                    self.diagonals.push_front(new_ivals.remove(&new_diag)
                        .unwrap_or_else(|| VisitedIntervalTree::new()));
                }

                self.k_min = new_kmin;
            }

            if new_kmax > self.k_max {
                let kmax_diff = self.k_max.diff(&new_kmax);
                eprintln!("k_max old: {:?}, k_max new: {:?}, diff: {:?}", self.k_max, new_kmax, kmax_diff);

                for i in 0..kmax_diff {
                    let new_diag = self.k_max + i as i64 + 1;
                    eprintln!("adding {:?}", new_diag);
                    self.diagonals.push_back(new_ivals.remove(&new_diag)
                        .unwrap_or_else(|| VisitedIntervalTree::new()));
                }

                self.k_max = new_kmax;
            }

            for (diag, itree) in DiagRange::closed(self.k_min, self.k_max).into_iter().zip(&mut self.diagonals) {
                if let Some(new_ivals_diag) = new_ivals.get(&diag) {
                    let to_merge = [itree, new_ivals_diag];
                    eprintln!("Merge after extend:");
                    let merged_ivals = VisitedIntervalTree::from_iter(merge_intervals(&to_merge));

                    *itree = merged_ivals;
                }
            }
        }

        // Check if we reached the end state
        let end_offset = O::new(seq.len());
        for end_node in graph.end_nodes() {
            let end_row = graph.node_ix_to_row(*end_node);
            let end_diag = (end_row, end_offset).diag();
            if let Some(ix) = self.diag_to_ix(end_diag) {
                if let Some(end_ival) = self.diagonals[ix].find(end_offset, O::max_value()).next() {
                    return Some((end_diag, end_ival));
                }
            }
        }

        eprintln!("After extend: {:?}", self);

        None
    }
}
