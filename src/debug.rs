use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;
use std::sync::mpsc;
use std::thread::JoinHandle;

use crate::errors::PoastaError;

pub mod messages {
    use serde::{Serialize, Deserialize};
    use petgraph::graph::IndexType;
    use crate::aligner::offsets::OffsetType;
    use crate::aligner::scoring::AlignmentStateTree;
    use crate::aligner::state::TreeIndexType;
    use crate::graphs::NodeIndexType;
    use crate::graphs::poa::POAGraph;

    #[derive(Debug, Serialize, Deserialize)]
    pub enum DebugOutputMessage {
        Empty,
        NewSequence { seq_name: String, seq_length: usize, max_rank: usize },
        IntermediateGraph { graph_dot: String },
        AlignStateTree { tree_gml: String },
        Terminate,
    }

    impl DebugOutputMessage {
        pub fn new_from_graph<GIx>(graph: &POAGraph<GIx>) -> Self
        where
            GIx: IndexType
        {
            Self::IntermediateGraph { graph_dot: format!("{}", graph) }
        }

        pub fn new_from_state_tree<T, N, O, Ix>(tree: &T) -> Self
        where
            T: AlignmentStateTree<N, O, Ix>,
            N: NodeIndexType,
            O: OffsetType,
            Ix: TreeIndexType
        {
            let gml = format!("{}", tree);
            Self::AlignStateTree { tree_gml: gml }
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

fn write_msg(writer: &mut impl Write, msg: &messages::DebugOutputMessage) {
    match serde_json::to_string(&msg) {
        Ok(json) => {
            if let Err(e) = writeln!(writer, "{}", json) {
                eprintln!("Error writing message to debug output!\n{}", e);
            }
        },
        Err(e) => eprintln!("Could not serialize debug data to JSON!\n{}", e)
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

            for msg in receiver {
                match msg {
                    messages::DebugOutputMessage::Empty => (),
                    messages::DebugOutputMessage::NewSequence { ref seq_name , seq_length: _, max_rank: _} => {
                        let mut output_file = File::create(output_path.join(format!("{seq_name}.txt")))
                            .map(BufWriter::new)?;
                        curr_seq_name = seq_name.clone();

                        write_msg(&mut output_file, &msg);
                    },
                    messages::DebugOutputMessage::IntermediateGraph { graph_dot } => {
                        let fname = output_path.join(format!("graph_for_{}.dot", &curr_seq_name));
                        let mut dot_file = File::create(fname)?;
                        write!(dot_file, "{}", graph_dot)?
                    },
                    messages::DebugOutputMessage::AlignStateTree { tree_gml } => {
                        let fname = output_path.join(format!("aln_state_tree_for_{}.gml", &curr_seq_name));
                        let mut gml_file = File::create(fname)?;
                        write!(gml_file, "{}", tree_gml)?
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