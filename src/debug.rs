use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;
use std::sync::mpsc;
use std::thread::JoinHandle;

use crate::errors::PoastaError;

pub mod messages {
    use serde::{Serialize, Deserialize};
    use petgraph::graph::IndexType;
    use serde::de::DeserializeOwned;
    use crate::aligner::visited::AlignState;
    use crate::aligner::layers::Layer;
    use crate::aligner::offsets::OffsetType;
    use crate::graphs::poa::POAGraph;
    use serde_json;

    #[derive(Debug, Serialize, Deserialize)]
    pub enum DebugOutputMessage {
        Empty,
        NewSequence { seq_name: String, seq: String, max_rank: usize },
        IntermediateGraph { graph_dot: String },
        NewScore(i64),
        VisitedIntervals { align_state: AlignState, intervals_json: String },
        Terminate,
    }

    impl DebugOutputMessage {
        pub fn new_from_graph<GIx>(graph: &POAGraph<GIx>) -> Self
        where
            GIx: IndexType + DeserializeOwned
        {
            Self::IntermediateGraph { graph_dot: format!("{}", graph) }
        }

        pub(crate) fn new_from_layer<O: OffsetType>(align_state: AlignState, layer: &Layer<O>) -> Self {
            let intervals_json = match serde_json::to_string(layer) {
                Err(e) => {
                    eprintln!("Could not convert layer to JSON! {}", e);
                    String::new()
                },
                Ok(v) => v,
            };

            Self::VisitedIntervals {
                align_state,
                intervals_json
            }
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
            let mut output_log = None;

            for msg in receiver {
                match msg {
                    messages::DebugOutputMessage::Empty => (),
                    messages::DebugOutputMessage::NewSequence { ref seq_name , ref seq, max_rank} => {
                        let mut new_file = File::create(output_path.join(format!("{seq_name}.txt")))
                            .map(BufWriter::new)?;
                        curr_seq_name = seq_name.clone();

                        match writeln!(&mut new_file, "# seq_name: {}\n# seq: {}\n# length: {}\n# num_nodes: {}", 
                                       curr_seq_name, seq, seq.len(), max_rank) {
                            Ok(_) => (),
                            Err(e) => eprintln!("Could not write sequence metadata for {:?}: {e}", curr_seq_name)
                        }

                        output_log = Some(new_file);
                    },
                    messages::DebugOutputMessage::IntermediateGraph { graph_dot } => {
                        let fname = output_path.join(format!("graph_for_{}.dot", &curr_seq_name));
                        let mut dot_file = File::create(fname)?;
                        write!(dot_file, "{}", graph_dot)?
                    },
                    messages::DebugOutputMessage::NewScore(score) => {
                        if let Some(ref mut ofile) = output_log {
                            match writeln!(ofile, "# score: {}", score) {
                                Ok(_) => (),
                                Err(e) => eprintln!("Could not write new score: {}", e)
                            }
                        }
                    },
                    messages::DebugOutputMessage::VisitedIntervals { align_state, intervals_json } => {
                        if let Some(ref mut ofile) = output_log {
                            match writeln!(ofile, "# state: {:?}", align_state) {
                                Ok(_) => (),
                                Err(e) => eprintln!("Could not write intervals ({:?}): {}", align_state, e)
                            }

                            match writeln!(ofile, "{}", intervals_json) {
                                Ok(_) => (),
                                Err(e) => eprintln!("Could not write intervals ({:?}): {}", align_state, e)
                            }
                        }
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