use xdrfile::{Frame, Result, Trajectory, XTCTrajectory};

fn main() -> Result<()> {

    let path = std::env::args()
        .skip(1)
        .nth(0)
        .expect("please provide one xtc trajectory path");

    let mut trj = XTCTrajectory::open_read(path)?;

    let num_atoms = trj.get_num_atoms()?;
    let mut frame = Frame::with_len(num_atoms);

    let mut n = 0;
    while trj.read(&mut frame).is_ok() {
        n += 1;
    }
    eprintln!("xdrreader: read {n} frames");

    Ok(())
}
