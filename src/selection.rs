use std::num::NonZeroU64;

// Invariant: The selection is only valid if the frame it reads them into is appropriately sized.
// It is assumed that the frame is correctly sized, i.e.,
//     len(frame.atoms) == len(IndexList) == sum(Map) == Until
// but it is also fine if the frame is too large for the selection. Stated differently,
//     len(frame.atoms) >= len(IndexList)
//     len(frame.atoms) >= sum(Mask)
//     len(frame.atoms) >= Until
// However, for the section of the frame that is not accounted for by the selection, the output is
// undefined. This does not mean it is unsafe, but they cannot be interpreted as valid positions.
// For Map a further invariant exists:
//     len(Mask) <= len(encoded_atoms)
#[derive(Debug, Default, Clone)]
pub enum AtomSelection {
    /// Include all atoms.
    #[default]
    All,
    /// A list of the indices of the positions to include in the selection.
    ///
    /// Invariant: The indices are _unique_ and _consecutive_.
    IndexList(Vec<u32>),
    /// A mask of the positions to include in the selection.
    ///
    /// If the value of the mask at an index `n` is `true`, the position at that same index `n` is
    /// included in the selection.
    Mask(Vec<bool>), // TODO: Bitmap optimization?
    /// Index of the last position to be included in the selection.
    ///
    /// This is an inclusive stop value, such that a value of 8 will mean that a total of 8 atoms
    /// are read into the frame.
    Until(u32),
}

impl AtomSelection {
    /// Determine whether some index `idx` is included in this [`AtomSelection`].
    ///
    /// Will return [`None`] once the index is beyond the scope of this `AtomSelection`.
    pub fn is_included(&self, idx: usize) -> Option<bool> {
        let idx = idx as u32;
        match self {
            AtomSelection::All => Some(true),
            AtomSelection::IndexList(indices) => {
                if indices.last()? <= &idx {
                    None
                } else {
                    Some(indices.contains(&idx)) // TODO: This may be a very bad thing.
                }
            }
            AtomSelection::Mask(mask) => mask.get(idx as usize).copied(),
            AtomSelection::Until(until) => {
                if &idx <= until {
                    Some(true)
                } else {
                    None
                }
            }
        }
    }
}

#[derive(Debug, Default, Clone)]
pub enum FrameSelection {
    /// Include all frames that are in a trajectory.
    #[default]
    All,
    /// Include frames that lie within a certain [`Range`].
    Range(Range),
    /// Include frames that match the indices in this list.
    ///
    /// Invariant: The indices in the FrameList are _unique_ and _consecutive_.
    FrameList(Vec<usize>),
}

impl FrameSelection {
    /// Determine whether some index `idx` is included in this [`FrameSelection`].
    ///
    /// Will return [`None`] once the index is beyond the scope of this `FrameSelection`.
    pub fn is_included(&self, idx: usize) -> Option<bool> {
        match self {
            FrameSelection::All => Some(true),
            FrameSelection::Range(range) => range.is_included(idx as u64),
            FrameSelection::FrameList(indices) => {
                if *indices.last()? <= idx {
                    None
                } else {
                    Some(indices.contains(&idx)) // TODO: This may be a very bad thing.
                }
            }
        }
    }
}

/// A selection of [`Frame`](super::Frame)s to be read from an [`XTCReader`](super::XTCReader).
///
/// The `start` of a [`Selection`] is always bounded, and is zero by default.
/// The `end` may be bounded or unbounded. In case the end is unbounded ([`None`]), a `Selection`
/// instructs the `XTCReader` to just read up to and including the last frame. If it is bounded
/// by [`Some`] value, the frames up to that index will be read.
/// The `step` describes the number of frames that passed in each stride.
/// The number of skipped `Frame`s is equal to `step` - 1.
/// For instance, given a `step` of four, one `Frame` is read and the following three are skipped.
///
/// # Note
///
/// An instance where `start` > `end` is a valid `Selection`, but it will not make much sense,
/// since the `Selection` will be understood to produce zero steps.
#[derive(Debug, Clone, Copy)]
pub struct Range {
    /// The `start` of a [`Selection`] is always bounded, and is zero by default.
    pub start: u64,
    /// The `end` may be bounded or unbounded.
    ///
    /// In case the end is unbounded ([`None`]), a `Selection` instructs the `XTCReader` to just
    /// read up to and including the last frame. If it is bounded by [`Some`] value, the frames up
    /// to that index will be read.
    pub end: Option<u64>,
    /// The `step` describes the number of frames that passed in each stride.
    ///
    /// The number of skipped `Frame`s is equal to `step` - 1.
    /// For instance, given a `step` of four, one `Frame` is read and the following three are skipped.
    pub step: NonZeroU64,
}

impl Range {
    pub fn new(start: Option<u64>, end: Option<u64>, step: Option<NonZeroU64>) -> Self {
        let mut sel = Self {
            end,
            ..Self::default()
        };
        if let Some(start) = start {
            sel.start = start;
        }
        if let Some(step) = step {
            sel.step = step;
        }
        sel
    }

    /// Apply the selection to an [`ExactSizeIterator`].
    fn apply<T>(&self, it: impl ExactSizeIterator<Item = T>) -> impl ExactSizeIterator<Item = T> {
        let end = self.end.unwrap_or(it.len() as u64) as usize;
        it.take(end)
            .skip(self.start as usize)
            .step_by(self.step.get() as usize)
    }

    fn is_included(&self, idx: u64) -> Option<bool> {
        if let Some(end) = self.end {
            // TODO: Better syntax with some range idk?
            if end < idx {
                return None;
            }
        };
        let in_range = self.start < idx;
        let in_step = (idx + self.start) % self.step == 0;
        Some(in_range && in_step)
    }
}

impl Default for Range {
    fn default() -> Self {
        Self {
            start: 0,
            end: None,
            step: NonZeroU64::new(1).unwrap(),
        }
    }
}
