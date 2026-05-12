use clap::ValueEnum;

pub mod lasagna;
pub mod poasta;

#[cfg(feature = "test")]
pub mod poasta_vs_spoa;

/// An enum indicating what kind of alignment to perform. Shared between
/// `poasta align` and `lasagna`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum AlignmentSpan {
    /// Perform global alignment
    Global,

    /// Perform semi-global alignment, i.e., globally align query but allow free gaps in the graph
    /// at the beginning and end
    SemiGlobal,

    /// Perform ends-free alignment, i.e., indels at the beginning or end on either the query or
    /// graph are free
    EndsFree,
}

/// Log verbosity levels exposed on the CLI.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl LogLevel {
    pub fn as_tracing_str(self) -> &'static str {
        match self {
            LogLevel::Error => "error",
            LogLevel::Warn => "warn",
            LogLevel::Info => "info",
            LogLevel::Debug => "debug",
            LogLevel::Trace => "trace",
        }
    }
}

/// Cost model selection.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum CostModelKind {
    /// Affine gap cost (single piecewise-linear gap penalty).
    Affine,

    /// Linear gap cost (no gap-open penalty).
    Linear,

    /// Two-piece affine gap cost.
    TwoPiece,
}
