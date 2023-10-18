use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::mpsc;
use std::thread::JoinHandle;

use crate::errors::PoastaError;

pub mod messages {
    use std::error::Error;
    use std::io::BufWriter;
    use serde::{Serialize, Deserialize};
    use petgraph::graph::IndexType;
    use crate::aligner::offsets::OffsetType;
    use crate::aligner::state::{Score, StateGraph};
    use crate::graphs::NodeIndexType;
    use crate::graphs::poa::POAGraph;

    #[derive(Debug, Serialize, Deserialize)]
    pub enum DebugOutputMessage {
        Empty,
        NewSequence { seq_name: String, sequence: String, max_rank: usize },
        IntermediateGraph { graph_dot: String },
        StateGraph { graph_tsv: Vec<u8>, score: usize },
        Terminate,
    }

    impl DebugOutputMessage {
        pub fn new_from_graph<GIx>(graph: &POAGraph<GIx>) -> Self
        where
            GIx: IndexType
        {
            Self::IntermediateGraph { graph_dot: format!("{}", graph) }
        }

        pub fn new_from_sg<SG, N, O>(state_graph: &SG, score: Score) -> Result<Self, Box<dyn Error>>
        where
            SG: StateGraph<N, O>,
            N: NodeIndexType,
            O: OffsetType,
        {
            let mut buf_writer = BufWriter::new(Vec::default());
            state_graph.write_tsv(&mut buf_writer)?;

            let graph_tsv = buf_writer.into_inner()?;

            Ok(Self::StateGraph { graph_tsv, score: score.into() } )
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

            let mut curr_seq_name = "none".to_string();
            let mut curr_seq = "".to_string();
            let mut curr_max_rank = 0;

            for msg in receiver {
                match msg {
                    messages::DebugOutputMessage::Empty => (),
                    messages::DebugOutputMessage::NewSequence { seq_name , sequence, max_rank} => {
                        curr_seq_name = seq_name;
                        curr_seq = sequence;
                        curr_max_rank = max_rank;
                    },
                    messages::DebugOutputMessage::IntermediateGraph { graph_dot } => {
                        let fname = output_path.join(format!("graph_for_{}.dot", &curr_seq_name));
                        let mut dot_file = File::create(fname)?;
                        write!(dot_file, "{}", graph_dot)?;
                        dot_file.flush()?;
                    },
                    messages::DebugOutputMessage::StateGraph { graph_tsv, score } => {
                        let fname = output_path.join(format!("state_graph_for_{}.score{score}.tsv", &curr_seq_name));
                        let mut gml_file = File::create(fname)?;
                        writeln!(gml_file, "# seq: {curr_seq} - max rank: {curr_max_rank}")?;
                        gml_file.write_all(&graph_tsv)?;
                        gml_file.flush()?;
                    },
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