use molly::XTCReader;

fn main() -> std::io::Result<()> {
    let path = std::env::args()
        .skip(1)
        .nth(0)
        .expect("please provide one xtc trajectory path");

    let file = std::fs::File::open(path)?;
    let mut reader = XTCReader::new(file);

    let frames = reader.read_all_frames()?;
    let n = frames.len();
    eprintln!("readall: read {n} frames");

    Ok(())
}
