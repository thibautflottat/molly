//! Filter an xtc trajectory, quickly.
//!
//! By Marieke Westendorp, 2024.
//! <ma3ke.cyber@gmail.com>
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::num::{NonZeroU64, ParseIntError};
use std::path::PathBuf;
use std::str::FromStr;

use clap::Parser;
use molly::buffer::{Buffer, UnBuffered};
use molly::reader::NBYTES_POSITIONS_PRELUDE;
use molly::selection::{AtomSelection, FrameSelection, Range};
use molly::{padding, read_positions, Frame, Header, XTCReader};

fn frame_selection_parser(selection: &str) -> Result<FrameSelection, ParseIntError> {
    let mut components = selection.split(':');
    let start = components.next().map(|s| s.parse()).transpose()?;
    let end = components.next().map(|s| s.parse()).transpose()?;
    let step = components
        .next()
        .map(|s| NonZeroU64::from_str(s))
        .transpose()?;
    Ok(FrameSelection::Range(Range::new(start, end, step)))
}

fn atom_selection_parser(selection: &str) -> Result<AtomSelection, ParseIntError> {
    let until: u32 = selection.parse()?;
    Ok(AtomSelection::Until(until))
}

// TODO: Consider making this one of several subcommands. This one could be called something like
// `molly filter ...`. Another would be `molly info` or `molly summary` or something.
/// Filter an xtc trajectory according to frame and atom selections.
#[derive(Parser)]
struct Args {
    /// Input path (xtc).
    input: PathBuf,

    /// Output path (xtc).
    output: PathBuf,

    /// Frame selection in the format `start:stop:step`. Each of these values optional.
    ///
    // TODO: Make these examples into unit tests for the frame_selection_parser and its atom counterpart.
    // TODO: Verify that I didn't make any mistakes in these examples, once everything is up and running.
    /// - `:100` will select the first 100 frames.
    ///
    /// - `3:14` will select the 4th up to and including the 14th frames, 11 frames in total.
    ///
    /// - `:100:2` will select every second frame from the the first 100 frames, 50 in total.
    #[arg(short, long, value_parser=frame_selection_parser)]
    frame_selection: Option<FrameSelection>,

    // TODO: Consider explaining why this seemingly silly limitation exists. It may be confusing
    // to just drop it here, but explaining it is also quite the ride.
    /// Atom selection single `stop` value.
    ///
    /// For each frame that is read, the compressed positions up to the provided index will be
    /// stored into the output file.
    ///
    // TODO: Verify that I didn't make any mistakes in these examples, once everything is up and running.
    /// - `1312` selects the first 1312 frames.
    ///
    /// Note that according to the xtc format, when the number of atoms in the frame is equal to
    /// or less than 9 (natoms <= 9), the atoms will be stored in an uncompressed manner.
    #[arg(short, long, value_parser=atom_selection_parser)]
    atom_selection: Option<AtomSelection>,

    // TODO: Add some {on, off, auto} enum?
    /// Use non-buffered reading mode. (Reading mode is buffered by default.)
    ///
    /// This may be faster under some circumstances, especially when the atom selection includes
    /// most of the atoms in a frame.
    #[arg(long = "unbuffered", default_value_t=true, action=clap::ArgAction::SetFalse)]
    is_buffered: bool,

    /// Write the trajectory in reverse.
    ///
    /// Selection functions the same regardless of whether this flag is set.
    #[arg(long)]
    reverse: bool,

    /// Print the time value for the selected frames to standard output.
    #[arg(long)]
    times: bool,

    /// Print the step number for the selected frames to standard output.
    ///
    /// If both `times` and `steps` are active, they will be separated by tabs and printed in that order.
    #[arg(long)]
    steps: bool,
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();

    let mut writer = BufWriter::new(std::fs::File::create(args.output)?);
    let file = std::fs::File::open(args.input)?;
    let mut reader = XTCReader::new(file);

    let frame_selection = args.frame_selection.unwrap_or_default();
    let atom_selection = args.atom_selection.unwrap_or_default();
    filter_frames(
        &mut reader,
        args.is_buffered,
        &mut writer,
        &frame_selection,
        &atom_selection,
        args.reverse,
        args.times,
        args.steps,
    )
}

fn filter_frames(
    reader: &mut XTCReader<File>,
    is_buffered: bool,
    writer: &mut BufWriter<File>,
    frame_selection: &FrameSelection,
    atom_selection: &AtomSelection,
    reversed: bool,
    times: bool,
    steps: bool,
) -> std::io::Result<()> {
    let mut scratch = Vec::new();
    let offsets = reader.determine_offsets(frame_selection.until())?;
    let enumerated_offsets: Vec<_> = {
        let enumerated = offsets.iter().enumerate();
        if reversed {
            enumerated.rev().collect()
        } else {
            enumerated.collect()
        }
    };
    let mut stdout = std::io::stdout();
    let mut frame = Frame::default();
    for (idx, &offset) in enumerated_offsets {
        match frame_selection.is_included(idx) {
            Some(true) => {}
            Some(false) => continue,
            None if !reversed => continue, // If we are reversed, we can't just stop early.
            None => break,
        }

        // Go to the start of this frame.
        reader.file.seek(SeekFrom::Start(offset))?;

        // Start of by reading the header.
        let header = reader.read_header()?;

        if times || steps {
            if times {
                write!(stdout, "{:.3}\t", header.time)?;
            }

            if steps {
                write!(stdout, "{}", header.step)?;
            }

            writeln!(stdout, "")?;
            continue;
        }

        // Now, we read the atoms.
        let natoms_frame = header.natoms; // The number of atoms specified for the frame.
        let nbytes = if natoms_frame <= 9 {
            // In this case, the positions are uncompressed. Each consists of three f32s, so we're
            // done pretty quickly.
            reader.read_smol_positions(natoms_frame, &mut frame, atom_selection)?
        } else {
            let nbytes = match is_buffered {
                false => read_positions::<UnBuffered, File>(
                    &mut reader.file,
                    natoms_frame,
                    &mut scratch,
                    &mut frame,
                    atom_selection,
                )?,
                true => read_positions::<Buffer, File>(
                    &mut reader.file,
                    natoms_frame,
                    &mut scratch,
                    &mut frame,
                    atom_selection,
                )?,
            };
            reader.step += 1;
            nbytes
        };

        // The number of atoms we are actually interested in for our output. Important to know
        // since it may be the case that more atoms are selected than are in the frame.
        let natoms = frame.natoms();
        // Reset to the start of the frame again, and skip the header.
        let offset_and_header = offset + Header::SIZE as u64;
        reader.file.seek(SeekFrom::Start(offset_and_header))?;

        // Redefine the header to reflect our changes.
        let header = Header {
            natoms,
            natoms_repeated: natoms,
            ..header
        };
        // And write it.
        writer.write_all(&header.to_be_bytes())?;

        if natoms <= 9 {
            // The number of positions is small. We encode the positions as uncompressed floats.
            for pos in &frame.positions {
                writer.write_all(&pos.to_be_bytes())?;
            }
        } else {
            // TODO: Consider 're-using' the scratch buffer!! It will contain (more than) the bytes we want to write out!
            // TODO: Invent some sort of SCRATCH mechanism here again.

            // Just copy over the precision, prelude, followed by the section of compressed bytes.
            let mut precision = [0; 4];
            reader.file.read_exact(&mut precision)?;
            writer.write_all(&precision)?;

            // Copy over the prelude, since that remains exactly the same.
            let mut prelude = [0; NBYTES_POSITIONS_PRELUDE];
            reader.file.read_exact(&mut prelude)?;
            writer.write_all(&prelude)?;

            let mut nbytes_old = [0; 4];
            reader.file.read_exact(&mut nbytes_old)?;
            // Check whether we totally messed up.
            let nbytes_old = u32::from_be_bytes(nbytes_old);
            assert!(
                nbytes <= nbytes_old as usize,
                "the new number of bytes ({nbytes}) must never be greater than the old number of bytes ({nbytes_old})"
            );

            // Write the new number of upcoming bytes.
            writer.write_all(&(nbytes as u32).to_be_bytes())?;
            // Note that we are dealing with xdr padding, here! (32-bit blocks.)
            let mut bytes = vec![0; nbytes + padding(nbytes)];
            reader.file.read_exact(&mut bytes[..nbytes])?;
            writer.write_all(&bytes)?;
        }
    }

    Ok(())
}
