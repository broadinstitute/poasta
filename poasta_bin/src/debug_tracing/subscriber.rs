use std::fmt;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::{Path, PathBuf};

use poasta::aligner::astar::AlignState;
use tracing::span;
use tracing_subscriber::Layer;

pub struct AlignStateLayer {
    output_dir: PathBuf,
}

impl AlignStateLayer {
    pub fn new(output_dir: &Path) -> Self {
        AlignStateLayer {
            output_dir: output_dir.to_owned(),
        }
    }
}

impl<S> Layer<S> for AlignStateLayer
where
    S: tracing::Subscriber,
    S: for<'lookup> tracing_subscriber::registry::LookupSpan<'lookup>, 
{
    fn on_new_span(&self, attrs: &span::Attributes<'_>, id: &span::Id, ctx: tracing_subscriber::layer::Context<'_, S>) {
        let span = ctx.span(id).expect("Span not found");
        let mut fname = SeqId::default();
        attrs.record(&mut fname);
        
        if let Some(seq_id) = fname.0 {
            let visited_fname = self.output_dir.join(format!("{}_visited.tsv", seq_id));
            let queued_fname = self.output_dir.join(format!("{}_queued.tsv", seq_id));
            
            let mut file_visited = File::create(&visited_fname)
                .map(BufWriter::new)
                .expect("Could not create file");
            let mut file_queued = File::create(&queued_fname)
                .map(BufWriter::new)
                .expect("Could not create file");
            
            match writeln!(&mut file_visited, "score\tnode\tdiag\toffset\tstate") {
                Ok(_) => (),
                Err(e) => {
                    eprintln!("Could not open visited file {:?}: {}", visited_fname, e);
                    return;
                }
            }
            match writeln!(&mut file_queued, "score\tnode\tdiag\toffset\tstate\tpriority") {
                Ok(_) => (),
                Err(e) => {
                    eprintln!("Could not open queued file {:?}: {}", queued_fname, e);
                    return;
                }
            }
            
            let mut extensions = span.extensions_mut();
            extensions.insert::<SpanOutputFiles>(SpanOutputFiles { visited: file_visited, queued: file_queued });
        }
    }
    
    fn on_event(&self, event: &tracing::Event<'_>, ctx: tracing_subscriber::layer::Context<'_, S>) {
        if event.metadata().target().ends_with("set_visited") {
            let mut visited_state = VisitedAlignmentState::default();
            event.record(&mut visited_state);
            
            let Some(span) = ctx.event_span(event) else {
                eprintln!("No span found for {:?}", visited_state);
                return;
            };
            let mut extensions = span.extensions_mut();
            
            if let Some(files) = extensions.get_mut::<SpanOutputFiles>() {
                let visited = &mut files.visited;
                let result = writeln!(
                    visited, 
                    "{}\t{}\t{}\t{}\t{}", 
                    visited_state.score.unwrap_or(usize::MAX), 
                    visited_state.node.unwrap_or(usize::MAX), 
                    visited_state.diag.unwrap_or(isize::MAX), 
                    visited_state.offset.unwrap_or(usize::MAX), 
                    visited_state.state.unwrap_or(format!("{:?}", AlignState::Match))
                );
                
                match result {
                    Ok(_) => (),
                    Err(e) => {
                        eprintln!("Could not write to visited file: {}", e);
                        eprintln!("state: {:?}", visited)
                    }
                }
            }
        } else if event.metadata().target().ends_with("extend") {
            let mut extended_state = ExtendedState::default();
            event.record(&mut extended_state);
            
            let Some(span) = ctx.event_span(event) else {
                eprintln!("No span found for {:?}", extended_state);
                return;
            };
            let mut extensions = span.extensions_mut();
            
            if let Some(files) = extensions.get_mut::<SpanOutputFiles>() {
                let visited = &mut files.visited;
                
                let Some(ostart) = &extended_state.offset else {
                    return;
                };
                
                let Some(oend) = &extended_state.extended_offset else {
                    return;
                };
                
                
                let default_str = format!("{:?}", AlignState::Match);
                let state_str = format!("{:?}", &extended_state.state.as_ref().unwrap_or(&default_str));
                
                for visited_offset in *ostart+1..=*oend {
                    let result = writeln!(
                        visited, 
                        "{}\t{}\t{}\t{}\t{}", 
                        extended_state.score.unwrap_or(usize::MAX), 
                        extended_state.node.unwrap_or(usize::MAX), 
                        extended_state.diag.unwrap_or(isize::MAX), 
                        visited_offset, 
                        &state_str
                    );
                    
                    match result {
                        Ok(_) => (),
                        Err(e) => {
                            eprintln!("Could not write to visited file: {}", e);
                            eprintln!("state: {:?}", visited)
                        }
                    }
                }
            }
            
        } else if event.metadata().target().ends_with("queue_item") {
            
        }
    }
}

#[derive(Debug, Default)]
struct SeqId(Option<String>);

impl tracing::field::Visit for SeqId {
    fn record_str(&mut self, field: &tracing::field::Field, value: &str) {
        if field.name() == "seq_id" {
            self.0 = Some(String::from(value))
        } else {
            self.record_debug(field, &value as &dyn fmt::Debug)
        }
    }
    
    fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn std::fmt::Debug) {
        eprintln!("Ignoring field {} = {:?}", field.name(), value);
    }
}

struct SpanOutputFiles {
    visited: BufWriter<File>,
    queued: BufWriter<File>,
}

#[derive(Debug, Default)]
struct VisitedAlignmentState {
    score: Option<usize>,
    node: Option<usize>,
    diag: Option<isize>,
    offset: Option<usize>,
    state: Option<String>,
}

impl tracing::field::Visit for VisitedAlignmentState {
    fn record_u64(&mut self, field: &tracing::field::Field, value: u64) {
        match field.name() {
            "score" => self.score = Some(value as usize),
            "node" => self.node = Some(value as usize),
            "offset" => self.offset = Some(value as usize),
            _ => self.record_debug(field, &value),
        }
    }
    
    fn record_i64(&mut self, field: &tracing::field::Field, value: i64) {
        match field.name() {
            "diag" => self.diag = Some(value as isize),
            _ => self.record_debug(field, &value),
        }
    }
    
    fn record_str(&mut self, field: &tracing::field::Field, value: &str) {
        match field.name() {
            "state" => self.state = Some(String::from(value)),
            "message" => (),
            _ => self.record_debug(field, &value),
        }
    }
    
    fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn std::fmt::Debug) {
        match field.name() {
            "state" => self.state = Some(format!("{:?}", value)),
            "message" => (),
            _ => eprintln!("Ignoring field {} = {:?}", field.name(), value)
        }
    }
}

#[derive(Debug, Default)]
struct ExtendedState {
    score: Option<usize>,
    node: Option<usize>,
    diag: Option<isize>,
    offset: Option<usize>,
    extended_offset: Option<usize>,
    state: Option<String>,
}

impl tracing::field::Visit for ExtendedState {
    fn record_u64(&mut self, field: &tracing::field::Field, value: u64) {
        match field.name() {
            "score" => self.score = Some(value as usize),
            "node" => self.node = Some(value as usize),
            "offset" => self.offset = Some(value as usize),
            "extended_offset" => self.extended_offset = Some(value as usize),
            _ => self.record_debug(field, &value),
        }
    }
    
    fn record_i64(&mut self, field: &tracing::field::Field, value: i64) {
        match field.name() {
            "diag" => self.diag = Some(value as isize),
            _ => self.record_debug(field, &value),
        }
    }
    
    fn record_str(&mut self, field: &tracing::field::Field, value: &str) {
        match field.name() {
            "state" => self.state = Some(String::from(value)),
            "message" => (),
            _ => self.record_debug(field, &value),
        }
    }
    
    fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn std::fmt::Debug) {
        match field.name() {
            "state" => self.state = Some(format!("{:?}", value)),
            "message" => (),
            _ => eprintln!("Ignoring field {} = {:?}", field.name(), value)
        }
    }
}