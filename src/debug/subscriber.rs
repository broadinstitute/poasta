//! A `tracing` layer that dumps the A* search, one file pair per aligned
//! sequence, for plotting and debugging.
//!
//! The layer watches for spans carrying a `seq_id` field. When one opens it
//! creates `<seq_id>_visited.tsv` and `<seq_id>_queued.tsv` in the output
//! directory, and every A* trace event fired inside that span (however deeply
//! nested) is appended to one of them:
//!
//! | event target                | file    | meaning                                 |
//! |-----------------------------|---------|-----------------------------------------|
//! | `poasta::astar::set_visited` | visited | a state popped from the queue and expanded |
//! | `poasta::astar::extend`      | visited | a state reached by greedy match extension  |
//! | `poasta::astar::queue_item`  | queued  | a state pushed onto the queue              |
//!
//! Both files are tab separated with a header row: `score`, `node` (the graph
//! node index), `offset` (query offset) and `state`; the queued file adds
//! `priority` (`score + heuristic`).

use std::fmt;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use tracing::span;
use tracing_subscriber::Layer;

use crate::debug::sanitise_filename;

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
    fn on_new_span(
        &self,
        attrs: &span::Attributes<'_>,
        id: &span::Id,
        ctx: tracing_subscriber::layer::Context<'_, S>,
    ) {
        let mut seq_id = SeqId::default();
        attrs.record(&mut seq_id);
        let Some(seq_id) = seq_id.0 else {
            return;
        };
        let Some(span) = ctx.span(id) else {
            return;
        };

        let stem = sanitise_filename(&seq_id);
        let open = |suffix: &str, header: &str| -> Option<BufWriter<File>> {
            let path = self.output_dir.join(format!("{stem}_{suffix}.tsv"));
            let mut file = match File::create(&path) {
                Ok(f) => BufWriter::new(f),
                Err(e) => {
                    eprintln!("Could not create {}: {e}", path.display());
                    return None;
                }
            };
            if let Err(e) = writeln!(file, "{header}") {
                eprintln!("Could not write {}: {e}", path.display());
                return None;
            }
            Some(file)
        };

        let (Some(visited), Some(queued)) = (
            open("visited", "score\tnode\toffset\tstate"),
            open("queued", "score\tnode\toffset\tstate\tpriority"),
        ) else {
            return;
        };

        span.extensions_mut()
            .insert::<SpanOutputFiles>(SpanOutputFiles { visited, queued });
    }

    fn on_event(&self, event: &tracing::Event<'_>, ctx: tracing_subscriber::layer::Context<'_, S>) {
        let target = event.metadata().target();
        let kind = if target.ends_with("set_visited") {
            EventKind::Visited
        } else if target.ends_with("extend") {
            EventKind::Extend
        } else if target.ends_with("queue_item") {
            EventKind::Queued
        } else {
            return;
        };

        // Events fire inside the engine's own spans; the files hang on the
        // enclosing `seq_id` span, so walk outwards to find them.
        let Some(span) = ctx
            .event_span(event)
            .and_then(|s| s.scope().find(|s| s.extensions().get::<SpanOutputFiles>().is_some()))
        else {
            return;
        };
        let mut extensions = span.extensions_mut();
        let Some(files) = extensions.get_mut::<SpanOutputFiles>() else {
            return;
        };

        let mut state = AlignmentState::default();
        event.record(&mut state);

        let result = match kind {
            EventKind::Visited => writeln!(
                files.visited,
                "{}\t{}\t{}\t{}",
                state.score, state.node, state.offset, state.state
            ),
            EventKind::Extend => {
                let mut r = Ok(());
                for offset in state.offset + 1..=state.extended_offset {
                    r = writeln!(
                        files.visited,
                        "{}\t{}\t{}\t{}",
                        state.score, state.node, offset, state.state
                    );
                }
                r
            }
            EventKind::Queued => writeln!(
                files.queued,
                "{}\t{}\t{}\t{}\t{}",
                state.score, state.node, state.offset, state.state, state.priority
            ),
        };
        if let Err(e) = result {
            eprintln!("Could not write A* debug output: {e}");
        }
    }
}

enum EventKind {
    Visited,
    Extend,
    Queued,
}

struct SpanOutputFiles {
    visited: BufWriter<File>,
    queued: BufWriter<File>,
}

#[derive(Debug, Default)]
struct SeqId(Option<String>);

impl tracing::field::Visit for SeqId {
    fn record_str(&mut self, field: &tracing::field::Field, value: &str) {
        if field.name() == "seq_id" {
            self.0 = Some(String::from(value))
        }
    }

    fn record_debug(&mut self, _field: &tracing::field::Field, _value: &dyn fmt::Debug) {}
}

/// Fields shared by all three event kinds; unused ones stay at their default.
#[derive(Debug, Default)]
struct AlignmentState {
    score: u64,
    node: u64,
    offset: u64,
    extended_offset: u64,
    priority: u64,
    state: String,
}

impl tracing::field::Visit for AlignmentState {
    fn record_u64(&mut self, field: &tracing::field::Field, value: u64) {
        match field.name() {
            "score" => self.score = value,
            "node" => self.node = value,
            "offset" => self.offset = value,
            "extended_offset" => self.extended_offset = value,
            "priority" => self.priority = value,
            _ => {}
        }
    }

    fn record_i64(&mut self, field: &tracing::field::Field, value: i64) {
        if value >= 0 {
            self.record_u64(field, value as u64)
        }
    }

    fn record_str(&mut self, field: &tracing::field::Field, value: &str) {
        if field.name() == "state" {
            self.state = value.to_owned();
        }
    }

    fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn fmt::Debug) {
        if field.name() == "state" {
            self.state = format!("{value:?}");
        }
    }
}
