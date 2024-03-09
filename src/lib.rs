#![feature(array_chunks, iter_array_chunks, array_try_map)]

use glam::{IVec3, Mat3, Vec3};

use crate::reader::{read_boxvec, read_compressed_positions, read_f32, read_f32s, read_i32};

mod reader;

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
    pub frames: Vec<Frame>,
    step: usize,
}

impl<R: std::io::Read> XTCReader<R> {
    pub const MAGIC: i32 = 1995;

    pub fn new(reader: R) -> Self {
        Self {
            file: reader,
            frames: Vec::new(),
            step: 0,
        }
    }

    /// Reads and returns a [`Frame`] and advances one step.
    pub fn read_frame(&mut self, frame: &mut Frame) -> std::io::Result<()> {
        let mut scratch = Vec::new();
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
        assert_eq!(magic, Self::MAGIC);
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
