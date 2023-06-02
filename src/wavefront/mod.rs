extern crate num;

use std::cmp;
use std::fmt::Debug;
use std::ops::{Index, IndexMut};

pub mod fr_points;
pub mod aligner;
pub mod compute;

use fr_points::{DiagIx, OffsetPrimitive, FRPoint, FRPointContainer, WFPoint};
use crate::wavefront::fr_points::{ExtendCandidate, OffsetWithBacktrace};


#[derive(Debug, Clone)]
struct Wavefront<Offset: OffsetPrimitive> {
    pub k_lo: DiagIx,
    pub k_hi: DiagIx,

    furthest_points: FRPointContainer<Offset>,
}

impl<Offset: OffsetPrimitive> Wavefront<Offset> {
    pub fn new() -> Self {
        Wavefront::default()
    }

    pub fn new_with_fr_points(k_lo: DiagIx, k_hi: DiagIx,
                              fr_points: FRPointContainer<Offset>) -> Self {
        // Check which diagonal has first/last furthest point, and use that as k_lo and k_hi
        // respectively
        let mut actual_k_lo = None;
        let mut actual_k_hi = None;

        for (k, offset) in (k_lo..=k_hi).zip(&fr_points) {
            if offset.is_some() {
                if actual_k_lo.is_none() {
                    actual_k_lo = Some(k);
                }

                actual_k_hi = Some(k);
            }
        }

        if let (Some(new_k_lo), Some(new_k_hi)) = (actual_k_lo, actual_k_hi) {
            let start = (new_k_lo - k_lo) as usize;
            let end = (new_k_hi - k_lo) as usize;

            Self {
                k_lo: new_k_lo,
                k_hi: new_k_hi,
                furthest_points: fr_points.into_iter()
                    .enumerate()
                    .filter_map(|(i, p)| if i >= start && i <= end { Some(p) } else { None })
                    .collect()
            }
        } else {
            Self::default()
        }
    }

    fn diag_to_ix(&self, diag: DiagIx) -> Option<usize> {
        eprintln!("to_ix, k_lo: {:?}, k_hi: {:?}, diag: {:?}, len: {:?}", self.k_lo, self.k_hi, diag, self.furthest_points.len());
        if diag >= self.k_lo && diag <= self.k_hi {
            Some((diag - self.k_lo) as usize)
        } else {
            None
        }
    }

    fn ensure_size(&mut self, diag: DiagIx) {
        let diff = if diag < self.k_lo {
            (self.k_lo - diag) as usize
        } else if diag > self.k_hi {
            (diag - self.k_hi) as usize
        } else {
            0
        };

        eprintln!("k_lo: {:?}, k_hi: {:?}, requested: {:?}, diff: {:?}", self.k_lo, self.k_hi, diag, diff);

        if diff != 0 {
            self.furthest_points.reserve(self.furthest_points.len() + diff);
            if diag < self.k_lo {
                for _ in 0..diff {
                    self.furthest_points.push_front(None);
                }

                self.k_lo -= diff as i64;
            } else {
                for _ in 0..diff {
                    self.furthest_points.push_back(None);
                }

                self.k_hi += diff as i64;
            }
        }
    }

    pub fn extend_candidate(&mut self, candidate: &ExtendCandidate<Offset>) {
        self.ensure_size(candidate.curr().diag());
        let ix = self.diag_to_ix(candidate.curr().diag()).unwrap();

        // Extend this candidate, increase the offset for this diagonal by one (done in
        // `make_offset_with_bt`)
        let offset_with_bt = candidate.make_offset_with_bt();

        match self.furthest_points[ix] {
            Some(ref mut v) => {
                if *v <= offset_with_bt {
                    v.offset = offset_with_bt.offset;

                    if offset_with_bt.prev.is_some() {
                        v.prev = offset_with_bt.prev;
                    }
                } else {
                    eprintln!("EXTEND {:?} - not overriding existing FR point {:?}", offset_with_bt, *v);
                }
            },
            None => self.furthest_points[ix] = Some(offset_with_bt)
        };
    }

    pub fn get(&self, k: DiagIx) -> Option<Offset> {
        if k >= self.k_lo && k <= self.k_hi {
            self.furthest_points[(k - self.k_lo) as usize].as_ref().map(|v| v.offset)
        } else {
            None
        }
    }

    pub fn iter(&self) -> impl Iterator<Item=FRPoint<Offset>> + '_ {
        self.furthest_points.iter()
            .zip(self.k_lo..=self.k_hi)
            .filter_map(|v| match v {
                (Some(offset), diag) => Some((diag, offset.offset)),
                _ => None
            })
    }
}

impl<Offset: OffsetPrimitive> Default for Wavefront<Offset> {
    fn default() -> Self {
        Wavefront {
            k_lo: 0,
            k_hi: 0,
            furthest_points: vec![None].into(),
        }
    }
}

impl<Offset: OffsetPrimitive> Index<DiagIx> for Wavefront<Offset> {
    type Output = Option<OffsetWithBacktrace<Offset>>;

    fn index(&self, index: DiagIx) -> &Self::Output {
        assert!(index >= self.k_lo && index <= self.k_hi);

        &self.furthest_points[(index - self.k_lo) as usize]
    }
}

impl<Offset: OffsetPrimitive> IndexMut<DiagIx> for Wavefront<Offset> {
    fn index_mut(&mut self, index: DiagIx) -> &mut Self::Output {
        assert!(index >= self.k_lo && index <= self.k_hi);

        &mut self.furthest_points[(index - self.k_lo) as usize]
    }
}



