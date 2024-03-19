#![feature(array_chunks, iter_array_chunks, array_try_map)]

use std::io::{self, SeekFrom};
use std::{cell::Cell, path::Path};

use glam::{Mat3, Vec3};

use crate::reader::{read_boxvec, read_compressed_positions, read_f32, read_f32s, read_i32};
use crate::selection::{AtomSelection, FrameSelection};

pub mod reader;
pub mod selection;

thread_local! {
    /// A scratch buffer to read encoded bytes into for subsequent decoding.
    static SCRATCH: Cell<Vec<u8>> = const { Cell::new(Vec::new()) };
}

pub type BoxVec = Mat3;

#[derive(Debug, Default, Clone)]
pub struct Frame {
    pub step: u32,
    /// Time in picoseconds.
    pub time: f32,
    pub boxvec: BoxVec,
    pub precision: f32,
    pub positions: Vec<f32>,
}

impl Frame {
    /// Returns an iterator over the coordinates stored in this [`Frame`].
    pub fn coords(&self) -> impl Iterator<Item = Vec3> + '_ {
        self.positions.array_chunks().map(|&c| Vec3::from_array(c))
    }
}

#[derive(Debug, Clone)]
pub struct XTCReader<R> {
    pub file: R,
    step: usize,
}

impl XTCReader<std::fs::File> {
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        Ok(Self::new(file))
    }
}

impl<R: io::Read> XTCReader<R> {
    pub const MAGIC: i32 = 1995;

    pub fn new(reader: R) -> Self {
        Self {
            file: reader,
            step: 0,
        }
    }

    /// A convenience function to read all frames in a trajectory.
    ///
    /// It is likely more efficient to use [`XTCReader::read_frame`] if you are only interested in
    /// the values of a single frame at a time.
    pub fn read_all_frames(&mut self) -> io::Result<Box<[Frame]>> {
        let mut frames = Vec::new();
        loop {
            let mut frame = Frame::default();
            if let Err(err) = self.read_frame(&mut frame) {
                match err.kind() {
                    // We have found the end of the file. No more frames, we're done.
                    io::ErrorKind::UnexpectedEof => break,
                    // Something else went wrong...
                    _ => Err(err)?,
                }
            }
            frames.push(frame);
        }
        Ok(frames.into_boxed_slice())
    }

    /// Reads and returns a [`Frame`] and advances one step.
    pub fn read_frame(&mut self, frame: &mut Frame) -> io::Result<()> {
        self.read_frame_with_selection(frame, &AtomSelection::All)
    }

    /// Reads and returns a [`Frame`] according to the [`AtomSelection`], and advances one step.
    pub fn read_frame_with_selection(
        &mut self,
        frame: &mut Frame,
        atom_selection: &AtomSelection,
    ) -> io::Result<()> {
        // Take the thread-local SCRATCH and use that while decoding the values.
        let mut scratch = SCRATCH.take();
        self.read_frame_with_scratch(frame, &mut scratch, atom_selection)
    }

    /// Reads and returns a [`Frame`] and advances one step, internally reading the compressed data
    /// into `scratch`.
    ///
    /// # Note
    ///
    /// This function performs the work of [`XTCReader::read_frame`], but leaves all allocations to
    /// the caller.
    ///
    /// The contents of `scratch` should not be depended upon! It just serves as a scratch buffer
    /// for the inner workings of decoding.
    ///
    /// In most cases, [`XTCReader::read_frame`] is more than sufficient. This function only serves
    /// to make specific optimization possible.
    pub fn read_frame_with_scratch(
        &mut self,
        frame: &mut Frame,
        scratch: &mut Vec<u8>,
        atom_selection: &AtomSelection,
    ) -> io::Result<()> {
        let file = &mut self.file;

        // Start of by reading the header.
        let magic = read_i32(file)?;
        assert_eq!(
            magic,
            Self::MAGIC,
            "found invalid magic number '{magic}' ({magic:#0x})"
        );
        let natoms: usize = read_i32(file)?
            .try_into()
            .map_err(|err| io::Error::other(format!("could not read natoms: {err}")))?;
        let step: u32 = read_i32(file)?
            .try_into()
            .map_err(|err| io::Error::other(format!("could not read step: {err}")))?;
        let time = read_f32(file)?;

        // Read the frame data.
        let boxvec = read_boxvec(file)?;
        let natoms_repeated = read_i32(file)?
            .try_into()
            .map_err(|err| io::Error::other(format!("could not read second natoms: {err}")))?;
        assert_eq!(natoms, natoms_repeated);

        // Now, we read the atoms.
        if natoms <= 9 {
            // In case the number of atoms is very small, just read their uncompressed positions.
            frame.positions.resize(natoms * 3, 0.0);
            let mut buf = [0.0; 9 * 3]; // We have at most 9 atoms, so we handle them on the stack.
            let buf = &mut buf[..natoms * 3];
            read_f32s(file, buf)?;
            frame.positions.truncate(0);
            frame.positions.extend(
                buf.array_chunks()
                    .enumerate()
                    .filter_map(|(idx, pos): (usize, &[f32; 3])| {
                        if atom_selection.is_included(idx).unwrap_or_default() {
                            Some(pos)
                        } else {
                            None
                        }
                    })
                    .flatten(),
            );
            // TODO: It is unclear to me what to do with the precision in this case. It is
            // basically invalid, or just irrelevant here, since we don't decode them. They were
            // never compressed to begin with.
        } else {
            // If the atom_selection specifies fewer atoms, we will only allocate up to that point.
            let natoms_selected = match atom_selection {
                AtomSelection::All => natoms,
                AtomSelection::Mask(mask) => mask.iter().filter(|&&include| include).count(),
                AtomSelection::Until(end) => *end as usize,
            };
            let natoms = usize::min(natoms, natoms_selected);

            frame.positions.resize(natoms * 3, 0.0);
            frame.precision = read_f32(file)?;
            read_compressed_positions(
                file,
                &mut frame.positions,
                frame.precision,
                scratch,
                atom_selection,
            )?;
            scratch.truncate(0); // FIXME: This is not actually necessary I think.
        }

        self.step += 1;

        frame.step = step;
        frame.time = time;
        frame.boxvec = boxvec;

        Ok(())
    }
}

impl<R: io::Read + io::Seek> XTCReader<R> {
    /// Reset the reader to its initial position.
    ///
    /// Go back to the first frame.
    pub fn home(&mut self) -> io::Result<()> {
        self.file.seek(SeekFrom::Start(0))?;
        self.step = 0;
        Ok(())
    }

    /// Returns the offsets from the headers in this [`XTCReader<R>`] from its current position.
    ///
    /// The last value points one byte after the last byte in the reader.
    ///
    /// If this function is called when the internal reader is not at its starting position, the
    /// frame offsets _from_ its position are determined. If you wish to determine the offsets from
    /// the initial reader position, call [`XTCReader::home`] before calling this function.
    ///
    /// # Errors
    ///
    /// This function will pass through any reader errors.
    pub fn determine_offsets_exclusive(&mut self, until: Option<usize>) -> io::Result<Box<[u64]>> {
        let file = &mut self.file;
        // Remember where we start so we can return to it later.
        let start_pos = file.stream_position()?;

        let mut offsets = Vec::new();

        while until.map_or(true, |until| offsets.len() < until) {
            match read_i32(file) {
                Ok(Self::MAGIC) => {}
                Ok(weird) => Err(io::Error::other(format!(
                    "found invalid magic number '{weird}' ({weird:#0x})"
                )))?,
                Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(err) => Err(err)?,
            };
            file.seek(SeekFrom::Current(84))?;
            let skip: u64 = read_i32(file)?
                .try_into()
                .map_err(|err| io::Error::other(format!("could not read frame size: {err}")))?;
            let padding = (4 - (skip as i64 % 4)) % 4; // FIXME: Why, and also, can we do this better?
            let offset = file.seek(SeekFrom::Current(skip as i64 + padding))?;
            offsets.push(offset);
        }

        // Return back to where we started.
        file.seek(SeekFrom::Start(start_pos))?;

        Ok(offsets.into_boxed_slice())
    }

    /// Returns the offsets of this [`XTCReader<R>`] from its current position.
    ///
    /// The last value points to the start of the last frame.
    ///
    /// If this function is called when the internal reader is not at its starting position, the
    /// frame offsets _from_ its position are determined. If you wish to determine the offsets from
    /// the initial reader position, call [`XTCReader::home`] before calling this function.
    ///
    /// # Errors
    ///
    /// This function will pass through any reader errors.
    pub fn determine_offsets(&mut self, until: Option<usize>) -> io::Result<Box<[u64]>> {
        let mut offsets = vec![0];
        let exclusive = self.determine_offsets_exclusive(until)?;
        offsets.extend(exclusive.iter().take(exclusive.len().saturating_sub(1)));
        Ok(offsets.into_boxed_slice())
    }

    /// Returns the frame sizes of this [`XTCReader<R>`].
    ///
    /// # Errors
    ///
    /// This function will pass through any reader errors.
    pub fn determine_frame_sizes(&mut self, until: Option<usize>) -> io::Result<Box<[u64]>> {
        let starts = self.determine_offsets_exclusive(until)?;
        let ends = starts.iter().clone().skip(1);
        Ok(starts
            .iter()
            .zip(ends)
            .map(|(s, e)| e - s)
            .collect::<Vec<_>>()
            .into_boxed_slice())
    }

    /// Seeks to offset, then reads and returns a [`Frame`] and advances one step.
    pub fn read_frame_at_offset(
        &mut self,
        frame: &mut Frame,
        offset: u64,
        atom_selection: &AtomSelection,
    ) -> io::Result<()> {
        self.file.seek(SeekFrom::Start(offset))?;
        self.read_frame_with_selection(frame, atom_selection)
    }

    /// Append [`Frame`]s to the `frames` buffer according to a [`Selection`].
    ///
    /// If successful, it will return the number of frames that were read.
    /// This can be useful since the selection itself is not enough to tell how many frames will
    /// actually be read.
    pub fn read_frames(
        &mut self,
        frames: &mut impl Extend<Frame>,
        frame_selection: &FrameSelection,
        atom_selection: &AtomSelection,
    ) -> io::Result<usize> {
        // TODO: Make this into a FrameSelection::until(&self) -> Option<usize>
        let until = match frame_selection {
            FrameSelection::All => None,
            FrameSelection::Range(range) => range.end.map(|end| end as usize),
            FrameSelection::FrameList(list) => Some(list.iter().max().copied().unwrap_or_default()),
        };
        let offsets = self.determine_offsets(until)?;
        let mut n = 0;
        for (idx, &offset) in offsets.iter().enumerate() {
            match frame_selection.is_included(idx) {
                Some(true) => {}
                Some(false) => continue,
                None => break,
            }
            let mut frame = Frame::default();
            self.read_frame_at_offset(&mut frame, offset, atom_selection)?;
            frames.extend(Some(frame));
            n += 1;
        }

        Ok(n)
    }
}
