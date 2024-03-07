use chemfiles::{Trajectory, Frame};

fn main() -> std::io::Result<()> {
    let path = std::env::args()
        .skip(1)
        .nth(0)
        .expect("please provide one xtc trajectory path");

    let mut trajectory = Trajectory::open(path, 'r').unwrap();
    let mut frame = Frame::new();

    let mut n = 0;
    while trajectory.read(&mut frame).is_ok() {n+=1}
    // eprintln!("read {n} frames");
    // eprintln!("{} atoms", frame.positions().len());

    Ok(())
}
