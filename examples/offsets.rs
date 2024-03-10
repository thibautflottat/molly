use molly::XTCReader;

fn main() -> std::io::Result<()> {
    let path = std::env::args()
        .skip(1)
        .nth(0)
        .expect("please provide one xtc trajectory path");

    let file = std::fs::File::open(path)?;
    let mut reader = XTCReader::new(file);

    let start = std::time::Instant::now();
    let offsets = reader.determine_offsets()?;
    let end = std::time::Instant::now();
    let duration = (end - start).as_secs_f32() * 1000.0;
    let n = offsets.len();
    eprintln!("offsets: found {n} offsets");
    eprintln!("         took {duration:.6} ms");

    for offset in offsets.iter() {
        println!("{offset}")
    }

    Ok(())
}
