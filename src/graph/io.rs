use std::str;
use std::io::Write;

use petgraph::graph::IndexType;
use petgraph::visit::{IntoEdgeReferences, EdgeRef};

use crate::aligner::astar::AlignableGraphRef;
use crate::errors::PoastaIOError;
use super::poa::POASeqGraph;


fn graphviz_node_color(label: u8) -> &'static str {
    match label {
        b'A' => "#80BC42",
        b'C' => "#006DB6",
        b'G' => "#F36C3E",
        b'T' => "#B12028",
        _ => "#939393",
    }
}


pub fn graph_to_dot<Ix>(writer: &mut impl Write, graph: &POASeqGraph<Ix>) -> Result<(), PoastaIOError>
where
    Ix: IndexType
{
    let seq_names_str = graph
        .get_sequences()
        .iter()
        .map(|v| format!("{}:{}", v.name(), v.start_node().index()))
        .collect::<Vec<String>>()
        .join("\t");

    writeln!(writer, "# seq:\t{seq_names_str}")?;

    writeln!(writer, "digraph {{")?;
    writeln!(writer, "rankdir=\"LR\"")?;
    writeln!(
        writer,
        "node [shape=none, style=filled, fillcolor=\"#e3e3e3\", penwidth=0]"
    )?;
    writeln!(writer)?;

    for n in graph.get_nodes() {
        let mut node_str = vec![
            format!("{} [label=<<TABLE border=\"0\" cellspacing=\"5\" cellpadding=\"5\"><TR>", n.index())
        ];
        
        for (i, base) in graph.node_seq(n).iter().enumerate() {
            let base_char = char::from(*base);
            let base_cell = format!("<TD port=\"{}\"><font color=\"{}\">{}</font></TD>", i, graphviz_node_color(*base), base_char);
            node_str.push(base_cell);
        }
        
        node_str.push(String::from("</TR></TABLE>>]"));
        
        writeln!(
            writer,
            "{}",
            node_str.join("\n")
        )?;
    }

    for n in graph.get_nodes() {
        for aligned_ival in graph.graph[n].aligned_intervals.iter() {
            writeln!(writer, "{}:{} -> {}:{} [style=dotted, color=\"#939393\", constraint=false]",
                n.index(),
                aligned_ival.node_start(),
                aligned_ival.other().index(),
                aligned_ival.other_start()
            )?;
            
            writeln!(writer, "{}:{} -> {}:{} [style=dotted, color=\"#939393\", constraint=false]",
                n.index(),
                aligned_ival.node_end(),
                aligned_ival.other().index(),
                aligned_ival.other_end()
            )?;
        }
    }

    let max_num_seq = graph
        .graph
        .edge_references()
        .map(|e| e.weight().sequence_ids.len())
        .max()
        .unwrap_or(1);
    let min_weight = 1.0;
    let max_weight = 40.0;
    let min_penwidth = 0.5;
    let max_penwidth = 3.5;

    for e in graph.graph.edge_references() {
        let seq_list_str = e
            .weight()
            .sequence_ids
            .iter()
            .map(|v| format!("s{v}"))
            .collect::<Vec<String>>()
            .join(" ");

        let num_seq = e.weight().sequence_ids.len();
        let (scaled_weight, scaled_penwidth) = if max_num_seq > 0 {
            let scaled_weight = (min_weight
                + ((num_seq as f64 / max_num_seq as f64) * (max_weight - min_weight)))
                .round() as i64;
            let scaled_penwidth =
                min_penwidth + ((num_seq as f64 / max_num_seq as f64) * (max_penwidth - min_penwidth));
            
            (scaled_weight, scaled_penwidth)
            
        } else {
            (1, 1.0)
        };

        writeln!(
            writer,
            "{} -> {} [weight={}; penwidth={}; label={}; class=\"{}\"]",
            e.source().index(),
            e.target().index(),
            scaled_weight,
            scaled_penwidth,
            num_seq,
            seq_list_str
        )?;
    }

    writeln!(writer, "}}")?;
    Ok(())
}