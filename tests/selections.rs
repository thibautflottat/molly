use std::num::NonZeroU64;

use molly::selection::{AtomSelection as AS, FrameSelection as FS, Range};

mod common;
use common::trajectories;

const PATH: &str = trajectories::SMOL;

/// Read frames according to some [`FrameSelection`] and [`AtomSelection`] and return the
/// number of frames that were read.
///
/// After reading these, the `frames` buffer will be cleared and the reader returns home, such
/// that the reader is in its orginal state again.
fn count_frames(
    reader: &mut molly::XTCReader<std::fs::File>,
    frames: &mut Vec<molly::Frame>,
    frame_selection: FS,
    atom_selection: AS,
) -> std::io::Result<usize> {
    reader.read_frames::<true>(frames, &frame_selection, &atom_selection)?;
    let nframes = frames.len();
    frames.clear();
    reader.home()?;
    Ok(nframes)
}

/// Read frames according to some [`AtomSelection`] and return the number of frames that were
/// read.
///
/// After reading these, the `frames` buffer will be cleared and the reader returns home, such
/// that the reader is in its orginal state again.
fn count_atoms(
    reader: &mut molly::XTCReader<std::fs::File>,
    atom_selection: AS,
) -> std::io::Result<usize> {
    let mut frame = molly::Frame::default();
    reader.read_frame_with_selection(&mut frame, &atom_selection)?;
    reader.home()?;
    Ok(frame.coords().count())
}

/// Read frames according to some [`FrameSelection`] and [`AtomSelection`] and assert that the
/// number of frames that were read is equal to `expected`.
macro_rules! assert_frames {
    ($frame_selection:expr, $atom_selection:expr => $expected:expr) => {{
        let mut reader = molly::XTCReader::open(&PATH)?;
        let mut frames = Vec::new();
        let reader = &mut reader;
        let frames = &mut frames;
        assert_eq!(
            count_frames(reader, frames, $frame_selection, $atom_selection)?,
            $expected
        );
        Ok(())
    }};
}

/// Read a frame according to some [`AtomSelection`] and assert that the number of atoms that
/// were read is equal to `expected`.
macro_rules! assert_atoms {
    ($atom_selection:expr => $expected:expr) => {{
        let mut reader = molly::XTCReader::open(&PATH)?;
        let reader = &mut reader;
        assert_eq!(count_atoms(reader, $atom_selection)?, $expected);
        Ok(())
    }};
}

mod frame_selection {
    use std::collections::BTreeSet;

    use super::*;

    const NFRAMES: usize = 1001;

    /// All frames.
    #[test]
    fn all_frames() -> std::io::Result<()> {
        assert_frames!(FS::All, AS::All => NFRAMES)
    }
    /// FrameSelection should be independent of the atoms that are selected.
    #[test]
    fn all_frames_with_atom_selection() -> std::io::Result<()> {
        assert_frames!(FS::All, AS::Until(100) => NFRAMES)
    }

    /// This Range should be identical to FS::All.
    #[test]
    fn range_all() -> std::io::Result<()> {
        assert_frames!(FS::Range(Range::new(None, None, None)), AS::All => NFRAMES)
    }
    /// Read half of the frames.
    #[test]
    fn range_half() -> std::io::Result<()> {
        assert_frames!(FS::Range(Range::new(None, None, NonZeroU64::new(2))), AS::All => 501)
    }
    /// Read with a huge step.
    #[test]
    fn range_huge_step() -> std::io::Result<()> {
        assert_frames!(FS::Range(Range::new(None, None, NonZeroU64::new(2000))), AS::All => 1)
    }
    /// Read the last couple of frames.
    #[test]
    fn range_last_frames() -> std::io::Result<()> {
        assert_frames!(FS::Range(Range::new(Some(NFRAMES as u64 - 20), None, None)), AS::All => 20)
    }
    /// Read the last couple of frames with a step.
    #[test]
    fn range_last_frames_step() -> std::io::Result<()> {
        assert_frames!(FS::Range(Range::new(Some(981), None, NonZeroU64::new(3))), AS::All => 7)
    }
    /// Read a clamped range.
    #[test]
    fn range_clamped() -> std::io::Result<()> {
        assert_frames!(FS::Range(Range::new(Some(500), Some(750), None)), AS::All => 250)
    }
    /// Read a clamped range with a step.
    #[test]
    fn range_clamped_step() -> std::io::Result<()> {
        assert_frames!(FS::Range(Range::new(Some(500), Some(750), NonZeroU64::new(5))), AS::All => 50)
    }

    /// Read according to a list of indices.
    #[test]
    fn indices() -> std::io::Result<()> {
        assert_frames!(FS::FrameList(BTreeSet::from_iter([0, 1, 500])), AS::All => 3)
    }
    /// Read the first frame.
    #[test]
    fn indices_first_frame() -> std::io::Result<()> {
        assert_frames!(FS::FrameList(BTreeSet::from_iter([0])), AS::All =>1)
    }
    /// Read a single frame at some index.
    #[test]
    fn indices_single_frame() -> std::io::Result<()> {
        assert_frames!(FS::FrameList(BTreeSet::from_iter([100])), AS::All => 1)
    }
    /// Read only the last index.
    #[test]
    fn indices_last_frame() -> std::io::Result<()> {
        assert_frames!(FS::FrameList(BTreeSet::from_iter([NFRAMES - 1])), AS::All => 1)
    }
    /// Read just past the last index.
    #[test]
    fn indices_after_last_frame() -> std::io::Result<()> {
        assert_frames!(FS::FrameList(BTreeSet::from_iter([NFRAMES])), AS::All => 0)
    }
    /// Read according to a list of indices with some beyond the last frame.
    #[test]
    fn indices_within_range_and_outside() -> std::io::Result<()> {
        assert_frames!(FS::FrameList(BTreeSet::from_iter([0, 1, 500, NFRAMES * 2])), AS::All => 3)
    }
    /// Read according to an empty list.
    #[test]
    fn indices_empty_list() -> std::io::Result<()> {
        assert_frames!(FS::FrameList(BTreeSet::new()), AS::All => 0)
    }
}

mod atom_selection {
    use super::*;

    const NATOMS: usize = 24316;

    /// All atoms.
    #[test]
    fn all_atoms() -> std::io::Result<()> {
        assert_atoms!(AS::All => NATOMS)
    }

    /// Read until index 0.
    #[test]
    fn until_zero() -> std::io::Result<()> {
        assert_atoms!(AS::Until(0) => 0)
    }
    /// Read until index 1 to read the first atom.
    #[test]
    fn until_first() -> std::io::Result<()> {
        assert_atoms!(AS::Until(1) => 1)
    }
    /// Read until half the atoms.
    #[test]
    fn until_half() -> std::io::Result<()> {
        assert_atoms!(AS::Until(NATOMS as u32 / 2) => NATOMS / 2)
    }
    /// Read all atoms with a perfectly set until.
    #[test]
    fn until_up_to_end() -> std::io::Result<()> {
        assert_atoms!(AS::Until(NATOMS as u32) => NATOMS )
    }
    /// Read until just beyond the number of atoms.
    #[test]
    fn until_just_beyond() -> std::io::Result<()> {
        assert_atoms!(AS::Until(NATOMS as u32 + 1) => NATOMS)
    }
    /// Read until far beyond the number of atoms.
    #[test]
    fn until_far_beyond() -> std::io::Result<()> {
        assert_atoms!(AS::Until(NATOMS as u32 + 1000) => NATOMS)
    }

    /// Read according to a list of indices.
    #[test]
    fn indices() -> std::io::Result<()> {
        assert_atoms!(AS::from_index_list(&[0, 1, 500]) => 3)
    }
    /// Read according to an empty list.
    #[test]
    fn indices_empty_list() -> std::io::Result<()> {
        assert_atoms!(AS::from_index_list(&[]) => 0)
    }
    /// Read the first atom.
    #[test]
    fn indices_first_atom() -> std::io::Result<()> {
        assert_atoms!(AS::from_index_list(&[0]) => 1)
    }
    /// Read a single atom at some index.
    #[test]
    fn indices_single_atom() -> std::io::Result<()> {
        assert_atoms!(AS::from_index_list(&[100]) => 1)
    }
    /// Read only the last index.
    #[test]
    fn indices_last_atom() -> std::io::Result<()> {
        assert_atoms!(AS::from_index_list(&[NATOMS as u32 - 1]) => 1)
    }
    /// Read just beyond the last index.
    #[test]
    fn indices_after_last_atom() -> std::io::Result<()> {
        assert_atoms!(AS::from_index_list(&[NATOMS as u32]) => 0)
    }
    /// Read far beyond the last index.
    #[test]
    fn indices_far_beyond_last_atom() -> std::io::Result<()> {
        assert_atoms!(AS::from_index_list(&[NATOMS as u32 + 1000]) => 0)
    }
    /// Read according to a list of indices with some beyond the last atom.
    #[test]
    fn indices_within_range_and_outside() -> std::io::Result<()> {
        assert_atoms!(AS::from_index_list(&[0, 1, 500, NATOMS as u32 + 1000]) => 3)
    }

    /// Read according to a mask.
    #[test]
    fn mask() -> std::io::Result<()> {
        assert_atoms!(AS::Mask(vec![true, false, false, true, false, true]) => 3)
    }
    /// Read according to an empty mask.
    #[test]
    fn mask_empty_list() -> std::io::Result<()> {
        assert_atoms!(AS::Mask(vec![]) => 0)
    }
    /// Read the first atom.
    #[test]
    fn mask_first_atom() -> std::io::Result<()> {
        assert_atoms!(AS::Mask(vec![true]) => 1)
    }
    /// Read a single atom at some index.
    #[test]
    fn mask_single_atom() -> std::io::Result<()> {
        assert_atoms!(AS::Mask([vec![false; 100], vec![true]].concat()) => 1)
    }
    /// Read only the last index.
    #[test]
    fn mask_last_atom() -> std::io::Result<()> {
        let n = NATOMS;
        let mut mask = vec![false; n];
        mask[n - 1] = true;
        assert_atoms!(AS::Mask(mask) => 1)
    }
    /// Read just beyond the last index.
    #[test]
    fn mask_after_last_atom() -> std::io::Result<()> {
        let n = NATOMS + 1;
        let mut mask = vec![false; n];
        mask[n - 1] = true;
        assert_atoms!(AS::Mask(mask) => 0)
    }
    /// Read far beyond the last index.
    #[test]
    fn mask_far_beyond_last_atom() -> std::io::Result<()> {
        let n = NATOMS + 1000;
        let mut mask = vec![false; n];
        mask[n - 1] = true;
        assert_atoms!(AS::Mask(mask) => 0)
    }
    /// Read according to a list of mask with some beyond the last atom.
    #[test]
    fn mask_within_range_and_outside() -> std::io::Result<()> {
        let n = NATOMS + 1000;
        let mut mask = vec![false; n];
        mask[0] = true;
        mask[1] = true;
        mask[500] = true;
        mask[n - 500] = true;
        mask[n - 1] = true;
        assert_atoms!(AS::Mask(mask) => 3)
    }
}
