use molly::selection::{AtomSelection, FrameSelection, Range};
use molly::XTCReader;

fn main() -> std::io::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("please provide one xtc trajectory path");
    let range = args.next().unwrap_or(String::from("::"));
    let is_buffered = args.next().map(|s| s == "buffered").unwrap_or_default();

    let file = std::fs::File::open(path)?;
    let mut reader = XTCReader::new(file);

    let range = parse_frame_selection(&range);
    let frame_selection = FrameSelection::Range(range);
    let atom_selection = AtomSelection::All;
    let mut frames = Void;
    let n = match is_buffered {
        true => reader.read_frames::<true>(&mut frames, &frame_selection, &atom_selection)?,
        false => reader.read_frames::<false>(&mut frames, &frame_selection, &atom_selection)?,
    };
    eprintln!("reader: read {n} frames");

    Ok(())
}

fn parse_frame_selection(s: &str) -> Range {
    let mut components = s.split(':');
    let start = components.next().map(|s| s.parse().ok()).flatten();
    let end = components.next().map(|s| s.parse().ok()).flatten();
    let step = components
        .next()
        .map(|s| {
            s.parse::<u64>()
                .map(|s| s.try_into().expect("step size must be greater than zero"))
                .ok()
        })
        .flatten();
    Range::new(start, end, step)
}

struct Void;

impl<A> Extend<A> for Void {
    fn extend<T: IntoIterator<Item = A>>(&mut self, _iter: T) {}
}
