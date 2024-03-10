#![feature(array_chunks, iter_array_chunks, array_try_map)]

use std::io::SeekFrom;
use std::{cell::Cell, path::Path};

use glam::{IVec3, Mat3, Vec3};

use crate::reader::{read_boxvec, read_compressed_positions, read_f32, read_f32s, read_i32};
pub use crate::selection::Selection;

pub mod reader;
mod selection;

thread_local! {
    /// A scratch buffer to read encoded bytes into for subsequent decoding.
    static SCRATCH: Cell<Vec<u8>> = Cell::new({
    eprintln!("hey there");
    Vec::new()});
}

pub type BoxVec = Mat3;

#[derive(Debug, Default, Clone)]
pub struct Frame {
    pub step: u32,
    /// Time in picoseconds.
    pub time: f32,
    pub boxvec: BoxVec,
    pub precision: f32,
    pub positions: Vec<i32>,
}

impl Frame {
    pub fn positions<'f>(&'f self) -> impl Iterator<Item = f32> + 'f {
        let inv_precision = self.precision.recip();
        self.positions
            .iter()
            .map(move |&v| v as f32 * inv_precision)
    }

    pub fn coords<'f>(&'f self) -> impl Iterator<Item = Vec3> + 'f {
        let inv_precision = self.precision.recip();
        self.positions
            .array_chunks()
            .map(move |&c| IVec3::from_array(c).as_vec3() * inv_precision)
    }
}

#[derive(Debug, Clone)]
pub struct XTCReader<R> {
    pub file: R,
    step: usize,
}

impl XTCReader<std::fs::File> {
    pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let file = std::fs::File::open(path)?;
        Ok(Self::new(file))
    }
}

impl<R: std::io::Read> XTCReader<R> {
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
    pub fn read_all_frames(&mut self) -> std::io::Result<Box<[Frame]>> {
        let mut frames = Vec::new();
        loop {
            let mut frame = Frame::default();
            if let Err(err) = self.read_frame(&mut frame) {
                match err.kind() {
                    // We have found the end of the file. No more frames, we're done.
                    std::io::ErrorKind::UnexpectedEof => break,
                    // Something else went wrong...
                    _ => Err(err)?,
                }
            }
            frames.push(frame);
        }
        Ok(frames.into_boxed_slice())
    }

    /// Reads and returns a [`Frame`] and advances one step.
    pub fn read_frame(&mut self, frame: &mut Frame) -> std::io::Result<()> {
        // Take the thread-local SCRATCH and use that while decoding the values.
        let mut scratch = SCRATCH.take();
        self.read_frame_with_scratch(frame, &mut scratch)
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
    ) -> std::io::Result<()> {
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
            .expect("natoms must be a positive integer");
        let step: u32 = read_i32(file)?
            .try_into()
            .expect("step must be a positive integer");
        let time = read_f32(file)?;

        // Read the frame data.
        let boxvec = read_boxvec(file)?;
        let natoms_repeated = read_i32(file)?.try_into().unwrap();
        assert_eq!(natoms, natoms_repeated);

        // Now, we read the atoms.
        frame.positions.resize(natoms * 3, 0);
        if natoms <= 9 {
            // In case the number of atoms is very small, just read their uncompressed positions.
            let mut buf = [0.0; 9 * 3]; // We have at most 9 atoms, so we handle them on the stack.
            let buf = &mut buf[..natoms * 9];
            read_f32s(file, buf)?;

            // Note that we do some cursed trickery here. In order to be able to write a i32
            // positions buffer, we need to convert these direct floats to the i32s. A bit sad, but
            // this case is extremely uncommon.
            frame.precision = 1000.0;
            let positions = buf.iter().map(|v| (v * frame.precision) as i32);

            frame.positions.truncate(0);
            frame.positions.extend(positions);
        } else {
            let precision = read_f32(file)?;
            frame.precision = precision;
            read_compressed_positions(file, &mut frame.positions, scratch)?;
            scratch.truncate(0);
        }

        self.step += 1;

        frame.step = step;
        frame.time = time;
        frame.boxvec = boxvec;

        Ok(())
    }
}

impl<R: std::io::Read + std::io::Seek> XTCReader<R> {
    /// Returns the offsets of this [`XTCReader<R>`].
    ///
    /// The last value points one byte after the last byte in the reader.
    ///
    /// # Panics
    ///
    /// Panics if the value at an offset is not the [magic number][`Self::MAGIC`]. If this is the
    /// case, something seriously weird has occurred.
    ///
    /// # Errors
    ///
    /// This function will pass through any reader errors.
    pub fn determine_offsets_exclusive(&mut self) -> std::io::Result<Box<[u64]>> {
        let file = &mut self.file;
        // Remember where we start so we can return to it later.
        let start_pos = file.stream_position()?;

        let mut offsets = Vec::new();

        loop {
            match read_i32(file) {
                Ok(Self::MAGIC) => {}
                Ok(weird) => panic!("found invalid magic number '{weird}' ({weird:#0x})"),
                Err(err) if err.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(err) => Err(err)?,
            };
            file.seek(SeekFrom::Current(84))?;
            let skip: u64 = read_i32(file)?.try_into().unwrap();
            let padding = (4 - (skip as i64 % 4)) % 4; // FIXME: Why, and also, can we do this better?
            let offset = file.seek(SeekFrom::Current(skip as i64 + padding)).unwrap();
            offsets.push(offset);
        }

        // Return back to where we started.
        file.seek(SeekFrom::Start(start_pos))?;

        Ok(offsets.into_boxed_slice())
    }

    /// Returns the offsets of this [`XTCReader<R>`].
    ///
    /// The last value points to the start of the last frame.
    ///
    /// # Panics
    ///
    /// Panics if the value at an offset is not the [magic number][`Self::MAGIC`]. If this is the
    /// case, something seriously weird has occurred.
    ///
    /// # Errors
    ///
    /// This function will pass through any reader errors.
    pub fn determine_offsets(&mut self) -> std::io::Result<Box<[u64]>> {
        let mut offsets = vec![0];
        let exclusive = self.determine_offsets_exclusive()?;
        offsets.extend(exclusive.iter().take(exclusive.len().saturating_sub(1)));
        Ok(offsets.into_boxed_slice())
    }

    /// Returns the frame sizes of this [`XTCReader<R>`].
    ///
    /// # Errors
    ///
    /// This function will pass through any reader errors.
    pub fn determine_frame_sizes(&mut self) -> std::io::Result<Box<[u64]>> {
        let starts = self.determine_offsets_exclusive()?;
        let ends = starts.iter().clone().skip(1);
        Ok(starts
            .iter()
            .zip(ends)
            .map(|(s, e)| e - s)
            .collect::<Vec<_>>()
            .into_boxed_slice())
    }

    /// Reads and returns a [`Frame`] and advances one step.
    pub fn read_frame_at_offset(&mut self, frame: &mut Frame, offset: u64) -> std::io::Result<()> {
        self.file.seek(SeekFrom::Start(offset))?;
        self.read_frame(frame)
    }

    /// Append [`Frame`]s to the `frames` buffer according to a [`Selection`].
    ///
    /// If successful, it will return the number of frames that were read.
    /// This can be useful since the selection itself is not enough to tell how many frames will
    /// actually be read.
    pub fn read_frames(
        &mut self,
        frames: &mut Vec<Frame>,
        selection: Selection,
    ) -> std::io::Result<usize> {
        let offsets = self.determine_offsets()?;
        let selected_offsets = selection.apply(offsets.iter());
        let n = selected_offsets.len();
        for &offset in selected_offsets {
            let mut frame = Frame::default();
            self.read_frame_at_offset(&mut frame, offset)?;
            frames.push(frame);
        }

        Ok(n)
    }
}
