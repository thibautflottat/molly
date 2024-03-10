use std::num::NonZeroU64;

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
pub struct Selection {
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

impl Selection {
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
    pub fn apply<T>(
        &self,
        it: impl ExactSizeIterator<Item = T>,
    ) -> impl ExactSizeIterator<Item = T> {
        let end = self.end.unwrap_or(it.len() as u64) as usize;
        it.take(end)
            .skip(self.start as usize)
            .step_by(self.step.get() as usize)
    }
}

impl Default for Selection {
    fn default() -> Self {
        Self {
            start: 0,
            end: None,
            step: NonZeroU64::new(1).unwrap(),
        }
    }
}
