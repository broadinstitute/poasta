//! Graph serialization to disk using serde

use std::io::{Read, Write};
use serde::de::DeserializeOwned;
use serde::Serialize;
use crate::errors::PoastaError;
use crate::graphs::poa::POAGraph;

use petgraph::graph::IndexType;

pub fn save_graph<Ix>(graph: &POAGraph<Ix>, out: impl Write) -> Result<(), PoastaError>
where
    Ix: IndexType + Serialize
{
    bincode::serialize_into(out, graph)?;

    Ok(())
}

pub fn load_graph<Ix>(reader: impl Read) -> Result<POAGraph<Ix>, PoastaError>
where
    Ix: IndexType + DeserializeOwned
{
    let graph: POAGraph<Ix> = bincode::deserialize_from(reader)?;

    Ok(graph)
}