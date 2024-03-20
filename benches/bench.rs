use std::io::ErrorKind;

use bencher::{benchmark_group, benchmark_main, Bencher};
use molly::{
    selection::{AtomSelection, FrameSelection},
    Frame, XTCReader,
};

benchmark_main!(reading);
benchmark_group!(reading, read_frame, read_frames);

const PATH: &str = "tests/trajectories/adk_oplsaa.xtc";

fn read_frame(b: &mut Bencher) {
    let mut reader = XTCReader::open(PATH).unwrap();
    let mut frame = Frame::default();
    b.iter(|| {
        match reader.read_frame(&mut frame) {
            Ok(_) => {}
            Err(err) if err.kind() == ErrorKind::UnexpectedEof => reader.home().unwrap(),
            Err(err) => panic!("{err}")
        }
    });
}

fn read_frames(b: &mut Bencher) {
    let mut reader = XTCReader::open(PATH).unwrap();
    let mut frames = Vec::new();
    b.iter(|| {
        reader
            .read_frames(&mut frames, &FrameSelection::All, &AtomSelection::All)
            .unwrap();
    });
}
