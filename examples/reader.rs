use molly::XTCReader;

fn main() -> std::io::Result<()> {
    let path = std::env::args()
        .skip(1)
        .nth(0)
        .expect("please provide one xtc trajectory path");

    let file = std::fs::File::open(path)?;
    let mut reader = XTCReader::new(file);
    let mut frame = molly::Frame::default();

    let mut n = 0;
    let mut scratch = Vec::new();
    while reader.read_frame_with_scratch(&mut frame,  &mut scratch).is_ok() {
        n += 1;
    }
    eprintln!("reader: read {n} frames");

    Ok(())
}
