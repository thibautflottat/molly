#![feature(array_chunks, iter_array_chunks, array_try_map)]

use glam::{Mat3, Vec3};

use crate::reader::{read_boxvec, read_compressed_floats, read_f32, read_f32s, read_i32};

mod reader;

pub type BoxVec = Mat3;

#[derive(Debug, Default, Clone)]
pub struct Frame {
    pub step: u32,
    /// Time in picoseconds.
    pub time: f32,
    pub boxvec: BoxVec,
    pub positions: Box<[Vec3]>,
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
    pub fn read_frame(&mut self) -> std::io::Result<Frame> {
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
        let positions = if natoms <= 9 {
            // In case the number of atoms is very small, just read their uncompressed positions.
            read_f32s(file, natoms * 3)?.collect()
        } else {
            let precision = read_f32(file)?;
            read_compressed_floats(file, natoms * 3, precision)?
        };
        let atoms = Vec::from_iter(positions.array_chunks().cloned().map(Vec3::from_array));

        self.step += 1;

        Ok(Frame {
            step,
            time,
            boxvec,
            positions: atoms.into_boxed_slice(),
        })
    }
}
