use std::fmt::Debug;

pub mod extend;
pub mod offsets;
pub mod aligner;
pub mod compute;

use offsets::{OffsetPrimitive, OffsetContainer, OffsetCell, Backtrace};


/// A graph wavefront holds the furthest reaching query offsets for the nodes
/// traversed so far at a given alignment score.
///
/// Certain alignment scores are unattainable with a given set of alignment costs. Those wavefronts
/// will never have any offsets stored. This is reflected in our enum, with a special Null variant.
///
/// If offsets are available for a given score, then they are stored in a `Vec` in
/// topological order, i.e., the offset at index 0 belongs to the node with rank 0.
#[derive(Clone, Debug)]
enum GraphWavefront<Offset: OffsetPrimitive> {
    Null,
    WithOffsets { offsets: OffsetContainer<Offset> }
}

impl<Offset: OffsetPrimitive> Default for GraphWavefront<Offset> {
    fn default() -> Self {
        Self::Null
    }
}

impl<Offset: OffsetPrimitive> GraphWavefront<Offset> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn initial() -> Self {
        Self::WithOffsets {
            offsets: vec![Some(OffsetCell::initial())]
        }
    }

    pub fn new_with_offsets(offsets: OffsetContainer<Offset>) -> Self {
        if offsets.iter().any(|v| v.is_some()) {
            Self::WithOffsets { offsets }
        } else {
            Self::Null
        }
    }

    pub fn resize(&mut self, max_rank: usize) {
        match self {
            Self::Null => (),
            Self::WithOffsets { offsets} => offsets.resize(max_rank + 1, None)
        }
    }

    pub fn get(&self, node_rank: usize) -> Option<Offset> {
        match self {
            Self::Null => None,
            Self::WithOffsets { offsets } => {
                if node_rank < offsets.len() {
                    offsets[node_rank].as_ref().map(|v| v.offset())
                } else {
                    None
                }
            }
        }
    }

    pub fn get_if_endpoint(&self, node_rank: usize) -> Option<Offset> {
        match self {
            Self::Null => None,
            Self::WithOffsets { offsets } => {
                if node_rank < offsets.len() {
                    // Only return offset if it's a current alignment "end point", to prevent
                    // expanding from nodes that were "extended"
                    offsets[node_rank].as_ref().and_then(|v| v.offset_if_endpoint())
                } else {
                    None
                }
            }
        }
    }

    pub fn get_backtrace(&self, node_rank: usize) -> Option<&Backtrace> {
        match self {
            Self::Null => None,
            Self::WithOffsets { offsets } => {
                if node_rank < offsets.len() {
                    offsets[node_rank].as_ref().map(|v| v.backtrace())
                } else {
                    None
                }
            }
        }
    }

    pub fn get_offsets(&self) -> Option<&OffsetContainer<Offset>> {
        match self {
            Self::Null => None,
            Self::WithOffsets { offsets } => Some(offsets)
        }
    }

    pub fn set(&mut self, node: usize, offset: OffsetCell<Offset>) {
        match self {
            Self::Null => panic!("Can't extend a point in a null-wavefront!"),
            Self::WithOffsets { offsets } => offsets[node] = Some(offset)
        }
    }

    pub fn iter_nodes(&self) -> WavefrontNodeIterator<'_, Offset> {
        WavefrontNodeIterator::new(self)
    }

    pub fn iter(&self) -> WavefrontOffsetIterator<'_, Offset> {
        WavefrontOffsetIterator::new(self)
    }

    pub fn iter_endpoints(&self) -> WavefrontEndpointIterator<'_, Offset> {
        WavefrontEndpointIterator::new(self)
    }


    pub fn len(&self) -> usize {
        match self {
            Self::Null => 0,
            Self::WithOffsets { offsets } => offsets.len()
        }
    }

    pub fn max_rank(&self) -> Option<usize> {
        match self {
            Self::Null => None,
            Self::WithOffsets { offsets } => Some(offsets.len() - 1)
        }
    }

}

/// Iterate over the node ranks that have an offset defined for a particular wavefront
struct WavefrontNodeIterator<'a, Offset: OffsetPrimitive> {
    /// Reference to the source wavefront
    wavefront: &'a GraphWavefront<Offset>,

    /// Current node rank
    curr: usize,

    /// To support iteration from the end we also have separate cursor starting from the back
    curr_back: usize
}

impl<'a, Offset: OffsetPrimitive> WavefrontNodeIterator<'a, Offset> {
    fn new(wavefront: &'a GraphWavefront<Offset>) -> Self {
        Self { wavefront, curr: 0, curr_back: wavefront.len() }
    }
}

impl<'a, Offset: OffsetPrimitive> Iterator for WavefrontNodeIterator<'a, Offset> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match self.wavefront {
            GraphWavefront::Null => None,
            GraphWavefront::WithOffsets { offsets } => {
                while self.curr < self.curr_back && offsets[self.curr].is_none() {
                    self.curr += 1;
                }

                if self.curr >= self.curr_back {
                    None
                } else {
                    let to_return = self.curr;
                    self.curr += 1;
                    Some(to_return)
                }
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.wavefront.len()))
    }
}

impl<'a, Offset: OffsetPrimitive> DoubleEndedIterator for WavefrontNodeIterator<'a, Offset> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match self.wavefront {
            GraphWavefront::Null => None,
            GraphWavefront::WithOffsets { offsets } => {
                if self.curr_back > 0 {
                    self.curr_back -= 1;
                }

                while self.curr_back > 0 && offsets[self.curr_back].is_none() {
                    self.curr_back -= 1;
                }

                if self.curr > 0 && self.curr_back < self.curr {
                    None
                } else {
                    // The while-loop above doesn't check if index 0 is Some
                    let to_return = offsets[self.curr_back].as_ref().map(|_| self.curr_back);

                    if self.curr_back > 0 {
                        self.curr_back -= 1;
                    } else {
                        // Trigger if-condition above to return None next iteration
                        self.curr += 1;
                    }

                    to_return
                }
            }
        }
    }
}


/// Iterate over (node rank, offset) tuples defined in a wavefront
struct WavefrontOffsetIterator<'a, Offset: OffsetPrimitive> {
    wavefront: &'a GraphWavefront<Offset>,

    /// Current node rank
    curr: usize,

    /// Separate cursor for double-ended iteration
    curr_back: usize
}

impl<'a, Offset: OffsetPrimitive> WavefrontOffsetIterator<'a, Offset> {
    fn new(wavefront: &'a GraphWavefront<Offset>) -> Self {
        Self { wavefront, curr: 0, curr_back: wavefront.len() }
    }
}

impl<'a, Offset: OffsetPrimitive> Iterator for WavefrontOffsetIterator<'a, Offset> {
    type Item = (usize, Offset);

    fn next(&mut self) -> Option<Self::Item> {
        match self.wavefront {
            GraphWavefront::Null => None,
            GraphWavefront::WithOffsets { offsets } => {
                while self.curr < self.curr_back && offsets[self.curr].is_none() {
                    self.curr += 1;
                }

                if self.curr >= self.curr_back {
                    None
                } else {
                    let to_return = offsets[self.curr].as_ref().map(|v| (self.curr, v.offset()));
                    self.curr += 1;

                    to_return
                }
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.wavefront.len()))
    }
}

impl<'a, Offset: OffsetPrimitive> DoubleEndedIterator for WavefrontOffsetIterator<'a, Offset> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match self.wavefront {
            GraphWavefront::Null => None,
            GraphWavefront::WithOffsets { offsets } => {
                if self.curr_back > 0 {
                    self.curr_back -= 1;
                }

                while self.curr_back > 0 && offsets[self.curr_back].is_none() {
                    self.curr_back -= 1;
                }

                if self.curr > 0 && self.curr_back < self.curr {
                    None
                } else {
                    // The while-loop above doesn't check if index 0 is Some
                    let to_return = offsets[self.curr_back].as_ref().map(|offset| (self.curr_back, offset.offset()));

                    if self.curr_back > 0 {
                        self.curr_back -= 1;
                    } else {
                        // Trigger if-condition above to return None next iteration
                        self.curr += 1;
                    }

                    to_return
                }
            }
        }
    }

}

/// Iterate over (node rank, offset) tuples defined in a wavefront, but only those that are an
/// alignment end-point
struct WavefrontEndpointIterator<'a, Offset: OffsetPrimitive> {
    wavefront: &'a GraphWavefront<Offset>,

    /// Current node rank
    curr: usize,

    /// Separate cursor for double-ended iteration
    curr_back: usize
}

impl<'a, Offset: OffsetPrimitive> WavefrontEndpointIterator<'a, Offset> {
    fn new(wavefront: &'a GraphWavefront<Offset>) -> Self {
        Self { wavefront, curr: 0, curr_back: wavefront.len() }
    }
}

impl<'a, Offset: OffsetPrimitive> Iterator for WavefrontEndpointIterator<'a, Offset> {
    type Item = (usize, Offset);

    fn next(&mut self) -> Option<Self::Item> {
        match self.wavefront {
            GraphWavefront::Null => None,
            GraphWavefront::WithOffsets { offsets } => {
                while self.curr < offsets.len() && offsets[self.curr].as_ref().map(|v| v.offset_if_endpoint()).is_none() {
                    self.curr += 1;
                }

                if self.curr >= self.curr_back {
                    None
                } else {
                    let offset = offsets[self.curr].as_ref().and_then(|v| v.offset_if_endpoint());
                    let to_return = offset.map(|v| (self.curr, v));
                    self.curr += 1;

                    to_return
                }
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.wavefront.len()))
    }
}

impl<'a, Offset: OffsetPrimitive> DoubleEndedIterator for WavefrontEndpointIterator<'a, Offset> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match self.wavefront {
            GraphWavefront::Null => None,
            GraphWavefront::WithOffsets { offsets } => {
                if self.curr_back > 0 {
                    self.curr_back -= 1;
                }

                while self.curr_back > 0 && offsets[self.curr_back].as_ref().and_then(|cell| cell.offset_if_endpoint()).is_none() {
                    self.curr_back -= 1;
                }

                if self.curr > 0 && self.curr_back < self.curr {
                    None
                } else {
                    let offset = offsets[self.curr_back].as_ref().and_then(|cell| cell.offset_if_endpoint());
                    let to_return = offset.map(|v| (self.curr_back, v));

                    if self.curr_back > 0 {
                        self.curr_back -= 1;
                    } else {
                        // Trigger if-condition above to return None next iteration
                        self.curr += 1;
                    }

                    to_return
                }
            }
        }
    }

}
