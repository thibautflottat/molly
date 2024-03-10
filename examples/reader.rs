use molly::{Selection, XTCReader};

fn main() -> std::io::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("please provide one xtc trajectory path");
    let range = args.next();

    let file = std::fs::File::open(path)?;
    let mut reader = XTCReader::new(file);
    let mut frame = molly::Frame::default();

    let mut n = 0;
    if let Some(range) = range {
        let selection = parse_selection(&range);
        dbg!(selection);
        let mut frames = Vec::new();
        n = reader.read_frames(&mut frames, selection)?;
    } else {
        while reader.read_frame(&mut frame).is_ok() {
            n += 1;
        }
    }
    eprintln!("reader: read {n} frames");

    Ok(())
}

fn parse_selection(s: &str) -> Selection {
    let mut components = s.split(':');
    let start = components.next().map(|s| s.parse().ok()).flatten();
    let end = components.next().map(|s| s.parse().ok()).flatten();
    let step = components.next().map(|s| {
        s.parse::<u64>()
            .map(|s| s.try_into().expect("step size must be greater than zero"))
            .ok()
    }).flatten();
    Selection::new(start, end, step)
}
