use std::io::{BufReader, ErrorKind};

use bencher::{benchmark_group, benchmark_main, Bencher};
use molly::{
    reader,
    selection::{AtomSelection, FrameSelection},
    Frame, XTCReader,
};

benchmark_main!(reading, decoding);
benchmark_group!(
    reading,
    read_frame,
    read_frames,
    read_frames_buffered,
    read_frames_few_atoms,
    read_frames_few_atoms_buffered
);
benchmark_group!(decoding, read_compressed_positions);

const PATH: &str = "tests/trajectories/adk_oplsaa.xtc";

fn read_frame(b: &mut Bencher) {
    let mut reader = XTCReader::open(PATH).unwrap();
    let mut frame = Frame::default();
    b.iter(|| match reader.read_frame(&mut frame) {
        Ok(_) => {}
        Err(err) if err.kind() == ErrorKind::UnexpectedEof => reader.home().unwrap(),
        Err(err) => panic!("{err}"),
    });
}

fn read_frames(b: &mut Bencher) {
    let mut reader = XTCReader::open(PATH).unwrap();
    let mut frames = Vec::new();
    b.iter(|| {
        reader
            .read_frames::<false>(&mut frames, &FrameSelection::All, &AtomSelection::All)
            .unwrap();
    });
}

fn read_frames_buffered(b: &mut Bencher) {
    let mut reader = XTCReader::open(PATH).unwrap();
    let mut frames = Vec::new();
    b.iter(|| {
        reader
            .read_frames::<true>(&mut frames, &FrameSelection::All, &AtomSelection::All)
            .unwrap();
    });
}

fn read_frames_few_atoms(b: &mut Bencher) {
    let mut reader = XTCReader::open(PATH).unwrap();
    let mut frames = Vec::new();
    b.iter(|| {
        reader
            .read_frames::<false>(&mut frames, &FrameSelection::All, &AtomSelection::Until(10))
            .unwrap();
    });
}

fn read_frames_few_atoms_buffered(b: &mut Bencher) {
    let mut reader = XTCReader::open(PATH).unwrap();
    let mut frames = Vec::new();
    b.iter(|| {
        reader
            .read_frames::<true>(&mut frames, &FrameSelection::All, &AtomSelection::Until(10))
            .unwrap();
    });
}

fn read_compressed_positions(b: &mut Bencher) {
    let natoms = 125;
    // A hand-tweaked test frame, derived from `delinyah_smaller.xtc`. Describes 125 positions.
    let bytes = include_bytes!("../tests/trajectories/delinyah_tiny.xtc");
    let start = 60; // We skip the header.
    let position_bytes = &bytes[start..];

    let mut positions = vec![0.0; natoms * 3];
    let mut scratch = Vec::new();
    let precision = 1000.0;
    b.iter(|| {
        let mut data = BufReader::new(position_bytes);
        reader::read_compressed_positions::<molly::buffer::UnBuffered, _>(
            &mut data,
            &mut positions,
            precision,
            &mut scratch,
            &AtomSelection::Until(natoms as u32),
        )
        .unwrap()
    });
}
