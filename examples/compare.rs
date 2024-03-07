use chemfiles::{Frame, Trajectory};
use molly::XTCReader;

fn round_to(v: f32, decimals: u32) -> f32 {
    let d = f32::powi(10.0, decimals as i32);
    (v * d).round() / d
}

fn main() -> std::io::Result<()> {
    let path = std::env::args()
        .skip(1)
        .nth(0)
        .expect("please provide one xtc trajectory path");

    let file = std::fs::File::open(&path)?;
    let mut reader = XTCReader::new(file);
    let mut trajectory = Trajectory::open(&path, 'r').unwrap();
    let mut cfframe = Frame::new();
    let mut n = 0;
    let mut natoms = 0;
    while let Ok(frame) = reader.read_frame() {
        trajectory.read(&mut cfframe).unwrap();

        for (a, &b) in frame.positions.iter().zip(cfframe.positions()) {
            let a = a.to_array().map(|v| round_to(v, 3));
            let b = b.map(|v| v as f32 * 0.1).map(|v| round_to(v, 3));
            // eprintln!("a, b = {a:?}\t\t{b:?}");
            assert_eq!(a, b);
        }

        natoms = frame.positions.len();
        n += 1;
    }
    eprintln!("read {n} frames");
    assert_eq!(natoms, cfframe.positions().len());
    eprintln!("{} atoms", cfframe.positions().len());

    Ok(())
}
