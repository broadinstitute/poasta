//! Graph serialization to disk using serde

use std::fmt;
use std::io::{Read, Write};
use petgraph::dot::Dot;
use petgraph::graph::IndexType;
use crate::errors::PoastaError;
use crate::graphs::poa::{POAGraph, POAGraphWithIx};

pub fn save_graph(graph: &POAGraphWithIx, out: impl Write) -> Result<(), PoastaError> {
    bincode::serialize_into(out, graph)?;

    Ok(())
}

pub fn load_graph(reader: impl Read) -> Result<POAGraphWithIx, PoastaError> {
    let graph: POAGraphWithIx = bincode::deserialize_from(reader)?;

    Ok(graph)
}

pub fn format_as_dot<Ix: IndexType>(writer: &mut impl fmt::Write, graph: &POAGraph<Ix>) -> fmt::Result {
    let transformed = graph.graph.map(
        |ix, data|
            format!("{:?} ({:?})", char::from(data.symbol), ix.index()),
        |_, data|
            format!("{}, {:?}", data.weight, data.sequence_ids)
    );

    let dot = Dot::new(&transformed);

    writeln!(writer, "{}", dot)?;

    Ok(())
}