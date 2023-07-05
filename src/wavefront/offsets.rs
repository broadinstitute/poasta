use std::cmp::{Ordering, Ord};
use std::fmt::Debug;

use num::{FromPrimitive, PrimInt};

pub type DiagIx = i64;

pub trait OffsetPrimitive: FromPrimitive + PrimInt
    + TryInto<usize> + Default + Copy + Debug { }

impl OffsetPrimitive for u8 { }
impl OffsetPrimitive for u16 { }
impl OffsetPrimitive for u32 { }
impl OffsetPrimitive for u64 { }

pub type OffsetContainer<Offset> = Vec<Option<OffsetCell<Offset>>>;

#[derive(PartialEq, Debug, Clone)]
pub enum BacktraceState {
    Start,
    End,
    Extend,
    Mismatch,
    InsOpen,
    InsExt,
    InsClose,
    DelOpen,
    DelExt,
    DelClose,

    // Two-piece (convex) alignment scoring
    InsOpen2,
    InsExt2,
    InsClose2,
    DelOpen2,
    DelExt2,
    DelClose2
}

impl AsRef<BacktraceState> for BacktraceState {
    fn as_ref(&self) -> &BacktraceState {
        self
    }
}

#[derive(Debug, Clone)]
pub struct Backtrace {
    prev_node: usize,
    state: BacktraceState
}

impl Backtrace {
    pub fn new(prev_node: usize, state: BacktraceState) -> Self {
        Backtrace { prev_node, state }
    }

    pub fn initial() -> Self {
        Backtrace { prev_node: 0, state: BacktraceState::Start }
    }

    pub fn prev_node(&self) -> usize {
        self.prev_node
    }

    pub fn state(&self) -> &BacktraceState {
        &self.state
    }

    pub fn align_state(&self) -> AlignState {
        self.state().into()
    }
}

pub enum AlignState {
    Start,
    End,
    Match,
    Mismatch,
    Insertion,
    Deletion
}

impl<T: AsRef<BacktraceState>> From<T> for AlignState {
    fn from(value: T) -> Self {
        match value.as_ref() {
            BacktraceState::Start => AlignState::Start,
            BacktraceState::End => AlignState::End,
            BacktraceState::Mismatch => AlignState::Mismatch,
            BacktraceState::Extend => AlignState::Match,
            BacktraceState::InsOpen | BacktraceState::InsOpen2 => AlignState::Match,
            BacktraceState::InsExt | BacktraceState::InsExt2 => AlignState::Insertion,
            BacktraceState::InsClose | BacktraceState::InsClose2 => AlignState::Insertion,
            BacktraceState::DelOpen | BacktraceState::DelOpen2 => AlignState::Match,
            BacktraceState::DelExt | BacktraceState::DelExt2 => AlignState::Deletion,
            BacktraceState::DelClose | BacktraceState::DelClose2 => AlignState::Deletion,
        }
    }
}


#[derive(Debug, Clone)]
pub enum CellType {
    EndPoint,
    Extended
}

#[derive(Debug, Clone)]
pub struct OffsetCell<Offset: OffsetPrimitive> {
    offset: Offset,
    cell_type: CellType,
    backtrace: Backtrace,
}

impl<Offset: OffsetPrimitive> OffsetCell<Offset> {
    pub fn new_endpoint(offset: Offset, backtrace: Backtrace) -> Self {
        Self { offset, cell_type: CellType::EndPoint, backtrace }
    }

    pub fn new_extended(offset: Offset, backtrace: Backtrace) -> Self {
        Self { offset, cell_type: CellType::Extended, backtrace }
    }

    pub fn initial() -> Self {
        Self::new_endpoint(Offset::zero(), Backtrace::initial())
    }

    pub fn set_cell_type(&mut self, cell_type: CellType) {
        self.cell_type = cell_type
    }

    pub fn is_endpoint(&self) -> bool {
        match self.cell_type {
            CellType::EndPoint => true,
            CellType::Extended => false
        }
    }

    pub fn is_extended(&self) -> bool {
        match self.cell_type {
            CellType::EndPoint => false,
            CellType::Extended => true
        }
    }

    pub fn offset(&self) -> Offset {
        self.offset
    }

    pub fn offset_if_endpoint(&self) -> Option<Offset> {
        match self.cell_type {
            CellType::Extended => None,
            CellType::EndPoint => Some(self.offset)
        }
    }

    pub fn backtrace(&self) -> &Backtrace {
        &self.backtrace
    }
}

#[derive(Clone, Debug)]
pub struct PrevCandidate<Offset: OffsetPrimitive>(Offset, Backtrace);

impl<Offset: OffsetPrimitive> PrevCandidate<Offset> {
    pub fn new(offset: Offset, bt_state: Backtrace) -> Self {
        Self(offset, bt_state)
    }

    pub fn offset(&self) -> Offset {
        self.0
    }

    pub fn backtrace(&self) -> &Backtrace {
        &self.1
    }

    pub fn into_offset_with_bt(self, new_offset: Offset) -> OffsetCell<Offset> {
        OffsetCell::new_endpoint(new_offset, self.1)
    }
}

impl<Offset: OffsetPrimitive> Ord for PrevCandidate<Offset> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl<Offset: OffsetPrimitive> PartialOrd for PrevCandidate<Offset> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<Offset: OffsetPrimitive> PartialEq for PrevCandidate<Offset> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<Offset: OffsetPrimitive> Eq for PrevCandidate<Offset> {

}