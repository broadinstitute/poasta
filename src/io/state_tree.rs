//! Output alignment state trees to a file
//!
//! This module provides functions to write the alignment state tree to various
//! graph file formats including GML. GML is also a file format supported by Python's
//! NetworkX, which could be used to easily analyze and plot the state tree from
//! Python.


use std::fmt::{self, Write, Debug};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::AlignmentStateTree;
use crate::aligner::state::{Backtrace, TreeIndexType};
use crate::graphs::NodeIndexType;

pub fn write_tree_gml<T, N, O, Ix>(writer: &mut impl Write, tree: &T) -> fmt::Result
where
    T: AlignmentStateTree<N, O, Ix>,
    N: Debug + NodeIndexType,
    O: Debug + OffsetType,
    Ix: TreeIndexType,
{
    writeln!(writer, "graph [")?;
    writeln!(writer, "  directed 1")?;

    // Write tree nodes
    for i in 0..tree.num_nodes() {
        let ix = Ix::new(i);
        writeln!(writer, "  node [")?;
        writeln!(writer, "    id {:?}", i)?;
        writeln!(writer, "    label {}", i)?;
        writeln!(writer, "    graph_node_ix {}", tree.get_node(ix).node().index())?;
        writeln!(writer, "    offset {:?}", tree.get_node(ix).offset())?;
        writeln!(writer, "    state \"{:?}\"", tree.get_node(ix).state())?;
        writeln!(writer, "  ]")?;
    }

    // Write edges
    for i in 0..tree.num_nodes() {
        let ix = Ix::new(i);
        if let Some(bt) = tree.get_node(ix).backtrace() {
            match bt {
                Backtrace::Step(parent) => {
                    writeln!(writer, "  edge [")?;
                    writeln!(writer, "    source {:?}", parent)?;
                    writeln!(writer, "    target {:?}", ix)?;
                    writeln!(writer, "    type \"single_step\"")?;
                    writeln!(writer, "    ext_matching_nodes \"_networkx_list_start\"")?;
                    writeln!(writer, "  ]")?;
                },
                Backtrace::ClosedIndel(parent) => {
                    writeln!(writer, "  edge [")?;
                    writeln!(writer, "    source {:?}", parent)?;
                    writeln!(writer, "    target {:?}", ix)?;
                    writeln!(writer, "    type \"closed_indel\"")?;
                    writeln!(writer, "    ext_matching_nodes \"_networkx_list_start\"")?;
                    writeln!(writer, "  ]")?;
                },
            }
        }
    }

    writeln!(writer, "]")?;

    Ok(())
}
