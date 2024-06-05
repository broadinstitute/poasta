use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::mpsc;
use std::thread::JoinHandle;

use crate::errors::PoastaError;

pub mod messages {
    use petgraph::graph::IndexType;
    #[cfg(feature = "serialize")]
    use serde::{Deserialize, Serialize};
    use crate::aligner::astar::AstarVisited;
    use crate::aligner::offsets::OffsetType;
    use crate::graphs::NodeIndexType;
    use crate::graphs::poa::POAGraph;

    #[cfg_attr(feature = "serialize", derive(Serialize, Deserialize))]
    #[derive(Debug)]
    pub enum DebugOutputMessage {
        Empty,
        NewSequence { seq_name: String, sequence: String, max_rank: usize },
        IntermediateGraph { graph_dot: String },
        AstarData { visited_tsv: String },
        Terminate,
    }

    impl DebugOutputMessage {
        pub fn new_from_graph<GIx>(graph: &POAGraph<GIx>) -> Self
        where
            GIx: IndexType
        {
            Self::IntermediateGraph { graph_dot: format!("{}", graph) }
        }

        pub fn new_from_astar_data<D, N, O>(data: &D) -> Self
            where D: AstarVisited<N, O>,
                  N: NodeIndexType,
                  O: OffsetType,
        {
            let mut tsv = String::new();
            data.write_tsv(&mut tsv).unwrap();
            Self::AstarData { visited_tsv: tsv }
        }
    }
}

pub struct DebugOutputWriter {
    transmitter: mpsc::Sender<messages::DebugOutputMessage>,
    worker: DebugOutputWorker
}

impl DebugOutputWriter {
    pub fn init<T: AsRef<Path> + Send>(debug_output_dir: T) -> Self {
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
            std::fs::create_dir_all(output_path.join("astar_iterations"))?;

            let mut curr_seq_name = "none".to_string();
            let mut curr_seq = "".to_string();
            let mut curr_max_rank = 0;
            let mut curr_astar_iter = 0usize;

            for msg in receiver {
                match msg {
                    messages::DebugOutputMessage::Empty => (),
                    messages::DebugOutputMessage::NewSequence { seq_name , sequence, max_rank} => {
                        curr_seq_name = seq_name;
                        curr_seq = sequence;
                        curr_max_rank = max_rank;
                        curr_astar_iter = 0;
                    },
                    messages::DebugOutputMessage::IntermediateGraph { graph_dot } => {
                        let fname = output_path.join(format!("graph_for_{}.dot", &curr_seq_name));
                        let mut dot_file = File::create(fname)?;
                        write!(dot_file, "{}", graph_dot)?;
                        dot_file.flush()?;
                    },
                    messages::DebugOutputMessage::AstarData { visited_tsv } => {
                        let fname = output_path.join(format!("astar_iterations/{curr_seq_name}.iter{curr_astar_iter}.tsv"));
                        let mut tsv_file = File::create(fname)?;
                        writeln!(tsv_file, "# seq_name: {curr_seq_name} - seq: {curr_seq} - max_rank: {curr_max_rank}")?;
                        tsv_file.write_all(visited_tsv.as_bytes())?;
                        tsv_file.flush()?;
                        curr_astar_iter += 1;
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
