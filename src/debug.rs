use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;
use std::sync::mpsc;
use std::thread::JoinHandle;

use petgraph::dot::Dot;

use crate::errors::PoastaError;

pub mod messages {
    use serde::{Serialize, Deserialize};
    use petgraph::Graph;
    use crate::graph::POAGraph;
    use crate::wavefront::offsets::OffsetPrimitive;
    use crate::wavefront::GraphWavefront;

    #[derive(Debug, Serialize, Deserialize)]
    pub enum DebugOutputMessage {
        Empty,
        NewSequence { seq_name: String },
        ExtendPath { start: (usize, usize), path: Vec<(usize, usize)> },
        CurrWavefront { score: i64, wf_type: String, node_offsets: Vec<(usize, usize, bool)> },
        IntermediateGraph { graph: Graph<String, String> },
        Terminate,
    }

    impl DebugOutputMessage {
        pub fn new_extended_path<T: TryInto<usize>>(start_node: usize, offset: T, path: Vec<(usize, usize)>) -> Self {
            if let Ok(offset_as_usize) = offset.try_into() {
                 Self::ExtendPath { start: (start_node, offset_as_usize), path }
            } else {
                panic!("Could not convert offset!")
            }
        }

        pub fn new_from_wavefront<T: OffsetPrimitive>(score: i64, wf_type: &str, wf: &GraphWavefront<T>) -> Self {
            let points: Option<Vec<(usize, usize, bool)>> = wf.get_offsets()
                .map(|v| {
                    v.iter()
                        .enumerate()
                        .filter_map(|(rank, cell)| {
                            cell.as_ref().map(|inner| {
                                if let Ok(offset) = inner.offset().try_into() {
                                    (rank, offset, inner.is_endpoint())
                                } else {
                                    panic!("Could not convert offset!")
                                }
                            })
                        })
                        .collect()
                });

            if let Some(v) = points {
                Self::CurrWavefront { score, wf_type: wf_type.to_string(), node_offsets: v }
            } else {
                Self::Empty
            }
        }

        pub fn new_from_graph(graph: &POAGraph) -> Self {
            let transformed = graph.graph.map(
                |ix, data|
                    format!("{:?} ({:?})", char::from(data.symbol), graph.get_node_rank(ix)),
                |_, data|
                    format!("{}, {:?}", data.weight, data.sequence_ids)
            );

            Self::IntermediateGraph { graph: transformed }
        }
    }
}

pub struct DebugOutputWriter {
    transmitter: mpsc::Sender<messages::DebugOutputMessage>,
    worker: DebugOutputWorker
}

impl DebugOutputWriter {
    pub fn new<T: AsRef<Path> + Send>(debug_output_dir: T) -> Self {
        let (tx, rx) = mpsc::channel();

        Self { transmitter: tx, worker: DebugOutputWorker::new(debug_output_dir, rx) }
    }

    pub fn log(&self, msg: messages::DebugOutputMessage) {
        if let Err(e) = self.transmitter.send(msg) {
            eprintln!("Could not log debug message!\nCause:{}", e)
        }
    }

    pub fn join(self) -> Result<(), PoastaError> {
        self.worker.join()
    }
}

struct DebugOutputWorker {
    thread: JoinHandle<Result<(), PoastaError>>,
}

impl DebugOutputWorker {
    fn new<T: AsRef<Path> + Send>(debug_output_dir: T, receiver: mpsc::Receiver<messages::DebugOutputMessage>) -> Self {
        let output_path = debug_output_dir.as_ref().to_path_buf();
        Self { thread: std::thread::spawn(move || {
            eprintln!("DEBUG: log output directory {:?}", output_path);
            std::fs::create_dir_all(&output_path)?;

            let mut output_file = File::create(output_path.join("startup.txt"))
                .map(BufWriter::new)?;
            let mut curr_seq_name = "none".to_string();

            for msg in receiver {
                eprintln!("Got debug message: {:?}", msg);

                match msg {
                    messages::DebugOutputMessage::Empty => (),
                    messages::DebugOutputMessage::NewSequence { seq_name } => {
                        output_file = File::create(output_path.join(format!("{seq_name}.txt")))
                            .map(BufWriter::new)?;

                        curr_seq_name = seq_name;
                    },
                    messages::DebugOutputMessage::ExtendPath { start: _, path: _ } |
                    messages::DebugOutputMessage::CurrWavefront { score: _, wf_type: _, node_offsets: _ } => {
                        match serde_json::to_string(&msg) {
                            Ok(json) => {
                                if let Err(e) = writeln!(output_file, "{}", json) {
                                    eprintln!("Error writing message to debug output!\n{}", e);
                                }
                            },
                            Err(e) => eprintln!("Could not serialize debug data to JSON!\n{}", e)
                        }
                    },
                    messages::DebugOutputMessage::IntermediateGraph { graph } => {
                        let fname = output_path.join(format!("graph_for_{}.dot", &curr_seq_name));
                        let dot = Dot::new(&graph);
                        let mut dot_file = File::create(fname)?;
                        write!(dot_file, "{:?}", dot)?

                    }
                    messages::DebugOutputMessage::Terminate => break
                }
            }

            Ok(())
        })}
    }

    fn join(self) -> Result<(), PoastaError> {
        self.thread.join().unwrap()
    }
}