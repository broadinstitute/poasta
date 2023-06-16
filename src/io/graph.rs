//! Graph serialization to disk using serde

use std::io::{Read, Write};
use crate::errors::PoastaError;
use crate::graph::POAGraph;

pub fn save_graph(graph: &POAGraph, out: impl Write) -> Result<(), PoastaError> {
    bincode::serialize_into(out, graph)?;

    Ok(())
}

pub fn load_graph(reader: impl Read) -> Result<POAGraph, PoastaError> {
    let graph: POAGraph = bincode::deserialize_from(reader)?;

    Ok(graph)
}