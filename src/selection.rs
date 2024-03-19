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
/// A selection of atoms.
#[derive(Debug, Default, Clone)]
pub enum AtomSelection {
    /// Include all atoms.
    #[default]
    All,
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
    /// Create a boolean mask from a list of indices.
    pub fn from_index_list(indices: &[u32]) -> Self {
        let max = match indices.iter().max() {
            Some(&max) => max as usize + 1,
            None => return Self::Mask(Vec::new()),
        };
        let mut mask = Vec::with_capacity(max);
        mask.resize(max, false);

        for &idx in indices {
            mask[idx as usize] = true;
        }

        Self::Mask(mask)
    }

    /// Determine whether some index `idx` is included in this [`AtomSelection`].
    ///
    /// Will return [`None`] once the index is beyond the scope of this `AtomSelection`.
    pub fn is_included(&self, idx: usize) -> Option<bool> {
        let idx = idx as u32;
        match self {
            AtomSelection::All => Some(true),
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

/// A selection of [`Frame`]s.
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
                if *indices.last()? < idx {
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
    /// to that index will be read. So, when `end` is bounded, it is an exclusive bound.
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

    fn is_included(&self, idx: u64) -> Option<bool> {
        if let Some(end) = self.end {
            // Determine whether `idx` is already beyond the defined range.
            if end <= idx {
                return None;
            }
        }
        let in_range = self.start <= idx;
        let in_step = self.step.get() == 1 || (idx + self.start) % self.step == 0;
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

#[cfg(test)]
mod tests {
    use super::*;

    mod frame {
        use std::num::NonZeroU64;

        use super::{FrameSelection, Range};

        #[test]
        fn zero_selection() {
            let list_empty = FrameSelection::FrameList(vec![]);
            let list_zero = FrameSelection::FrameList(vec![0]);
            let range_empty = FrameSelection::Range(Range::new(None, Some(0), None));

            for idx in 0..1000 {
                assert!(list_empty.is_included(idx).is_none());
                if idx > 0 {
                    assert!(list_zero.is_included(idx).is_none());
                }
                assert!(range_empty.is_included(idx).is_none());
            }
        }

        #[test]
        fn first_n() {
            let n = 100;
            let step = NonZeroU64::new(17).unwrap();

            let list = FrameSelection::FrameList((0..=n).collect());
            let until = FrameSelection::Range(Range::new(None, Some(n as u64), None));
            let from_n = FrameSelection::Range(Range::new(Some(n as u64), None, None));
            let until_stepped = FrameSelection::Range(Range::new(None, Some(n as u64), Some(step)));
            let from_n_stepped =
                FrameSelection::Range(Range::new(Some(n as u64), None, Some(step)));
            let all = FrameSelection::All;

            for idx in 0..2 * n {
                if idx < n {
                    assert_eq!(list.is_included(idx), Some(true));
                    assert_eq!(until.is_included(idx), Some(true));
                    assert_eq!(
                        until_stepped.is_included(idx),
                        Some(idx as u64 % step.get() == 0),
                    );
                } else {
                    assert!(list.is_included(idx).is_none());
                    assert!(until.is_included(idx).is_none());
                    assert!(until_stepped.is_included(idx).is_none());
                }
                let from_n_included = idx >= n;
                assert_eq!(from_n.is_included(idx), Some(from_n_included));
                assert_eq!(
                    from_n_stepped.is_included(idx),
                    Some(from_n_included && (n as u64 + idx as u64) % step.get() == 0),
                );
                assert_eq!(all.is_included(idx), Some(true));
            }
        }
    }

    mod atom {
        use super::AtomSelection;

        #[test]
        fn zero_selection() {
            let m = 100;

            let mask_empty = AtomSelection::Mask(vec![]);
            let mask_false = AtomSelection::Mask(vec![false; m]);
            let list_empty = AtomSelection::from_index_list(&vec![]);
            let list_zero = AtomSelection::from_index_list(&vec![0]);
            let until_zero = AtomSelection::Until(0);

            for idx in 0..1000 {
                assert!(mask_empty.is_included(idx).is_none());
                if idx < m {
                    assert_eq!(mask_false.is_included(idx), Some(false));
                } else {
                    assert!(mask_false.is_included(idx).is_none());
                }
                assert!(list_empty.is_included(idx).is_none());
                if idx > 0 {
                    assert!(until_zero.is_included(idx).is_none());
                    assert!(list_zero.is_included(idx).is_none());
                } else {
                    assert_eq!(until_zero.is_included(idx), Some(true));
                    assert_eq!(list_zero.is_included(idx), Some(true));
                }
            }
        }

        #[test]
        fn first_n() {
            let n = 100;
            let mask = AtomSelection::Mask(vec![true; n]);
            let mask_trailing_false = AtomSelection::Mask([vec![true; n], vec![false; n]].concat());
            let list = AtomSelection::from_index_list(&(0..n as u32).collect::<Vec<_>>());
            let until = AtomSelection::Until(n as u32 - 1);
            let all = AtomSelection::All;

            for idx in 0..2 * n {
                if idx < n {
                    assert_eq!(mask.is_included(idx), Some(true));
                    assert_eq!(list.is_included(idx), Some(true));
                    assert_eq!(until.is_included(idx), Some(true));
                } else {
                    assert!(mask.is_included(idx).is_none());
                    assert!(list.is_included(idx).is_none());
                    assert!(until.is_included(idx).is_none());
                }
                assert_eq!(mask_trailing_false.is_included(idx), Some(idx < n));
                assert_eq!(all.is_included(idx), Some(true));
            }
        }
    }
}
