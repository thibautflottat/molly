use molly::XTCReader;

fn main() -> std::io::Result<()> {
    let path = std::env::args()
        .skip(1)
        .nth(0)
        .expect("please provide one xtc trajectory path");

    let file = std::fs::File::open(path)?;
    let mut reader = XTCReader::new(file);

    let mut n = 0;
    while reader.read_frame().is_ok() {
        n += 1;
    }
    // eprintln!("n = {n}");

    Ok(())
}
