//! Implementation of the classical PO-PO alignment algorithm, aligning
//! one partial order graph to another partial order graph.

use std::slice;

use crate::align::cost_models::affine::Affine;
use crate::align::cost_models::AlignmentCostModel;
use crate::align::traits::AlignmentEngine;
use crate::errors::PoastaError;
use crate::graph::poa::{POAGraph, IndexType};
use crate::graph::traits::{GraphBase, GraphWithNodeOrdering};


pub struct POPOAEngine {
    cost_model: Affine
}


impl POPOAEngine {

}


impl<Ix> AlignmentEngine<POAGraph<Ix>> for POPOAEngine
where
    Ix: IndexType
{
    type Graph = POAGraph<Ix>;
    type Success = ();
    type Error = PoastaError<Ix>;

    fn align(&self, graph: &Self::Graph, to_align: &POAGraph<Ix>) -> Result<Self::Success, Self::Error> {
        // Subtract one to omit end node
        let matrix_width = to_align.node_count() - 1;
        let matrix_height = graph.node_count() - 1;
        let mut matrix_values = vec![0; 3 * matrix_width * matrix_height].into_boxed_slice();
        
        // SAFETY: non-overlapping slices to the same contiguous array allocated above.
        let matrix_len = matrix_width * matrix_height;
        let mut I = unsafe { slice::from_raw_parts_mut(&mut matrix_values[0], matrix_len) };
        let mut D = unsafe { slice::from_raw_parts_mut(&mut matrix_values[matrix_len], matrix_len) };
        let mut M = unsafe { slice::from_raw_parts_mut(&mut matrix_values[2 * matrix_len], matrix_len) };
        
        // Initialize first col of I
        for rank in 1..matrix_height {
            let ix = matrix_width * rank;
            let (pred_rows, curr_row) = I.split_at_mut(ix);
            let node = graph.rank_to_node(rank);
            let min_pred_cost = graph.predecessors(node)
                .map(|pred| {
                    let pred_ix = matrix_width * graph.node_rank(pred);
                    pred_rows[pred_ix]
                })
                .min()
                .unwrap_or(0usize);
            
            let cost = if min_pred_cost == 0 {
                self.cost_model.gap_open() as usize + self.cost_model.gap_extend() as usize
            } else {
                min_pred_cost + self.cost_model.gap_extend() as usize
            };
            
            curr_row[0] = cost;
        }
        
        // Initialize first row of D
        for rank in 1..matrix_width {
            let (pred_cols, curr_col) = D.split_at_mut(rank);
            let node = to_align.rank_to_node(rank);
            let min_pred_cost = to_align.predecessors(node)
                .map(|pred| {
                    let pred_col = to_align.node_rank(pred);
                    pred_cols[pred_col]
                })
                .min()
                .unwrap_or(0usize);
            
            let cost = if min_pred_cost == 0 {
                self.cost_model.gap_open() as usize + self.cost_model.gap_extend() as usize
            } else {
                min_pred_cost + self.cost_model.gap_extend() as usize
            };
            
            curr_col[0] = cost;
        }
        
        // Initialize first row of M
        for rank in 1..matrix_width {
            M[rank] = D[rank]
        }
        
        // Initialize first col of M
        for rank in 1..matrix_height {
            let ix = rank * matrix_width;
            M[ix] = I[ix];
        }
        
        // Fill in rest of the matrices using the DP recursion
        for rank_row in 1..matrix_height {
            let node_row = graph.rank_to_node(rank_row);
            let row_ix = rank_row * matrix_width;
            
            // Split slices to satisfy borrow checker (only mutating curr row)
            let (i_pred_rows, i_curr_row) = I.split_at_mut(row_ix);
            let (d_pred_rows, d_curr_row) = D.split_at_mut(row_ix);
            let (m_pred_rows, m_curr_row) = M.split_at_mut(row_ix);
            
            for rank_col in 1..matrix_width {
                let node_col = to_align.rank_to_node(rank_col);
                let col_ix = row_ix + rank_col;
                
                let (i_pred_cols, i_curr_col) = i_curr_row.split_at_mut(rank_col);
                let (d_pred_cols, d_curr_col) = d_curr_row.split_at_mut(rank_col);
                let (m_pred_cols, m_curr_col) = m_curr_row.split_at_mut(rank_col);
                
                // Update I
                // Check for gap extension
                let min_pred_i = to_align.predecessors(node_col)
                    .map(|pred| {
                        let pred_col = to_align.node_rank(pred);
                        i_pred_cols[pred_col]
                    })
                    .min();
                // Check for gap open
                let max_pred_m = to_align.predecessors(node_col)
                    .map(|pred| {
                        let pred_col = to_align.node_rank(pred);
                        
                    });
            }
        }
        
    }
}
