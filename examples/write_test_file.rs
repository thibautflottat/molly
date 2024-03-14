use xdrfile::{Frame, Trajectory, XTCTrajectory};

fn main() -> std::io::Result<()> {
    let mut args = std::env::args().skip(1);

    let path = args.next().expect("path");
    let natoms: usize = args.next().expect("natoms").parse().unwrap();
    let nframes: usize = args.next().expect("nframs").parse().unwrap();

    let mut traj = XTCTrajectory::open_write(path).unwrap();

    let frames = (0..nframes).map(|fi| {
        let mut frame = Frame::new();
        frame.step = fi;
        frame.time = fi as f32;
        frame.box_vector = [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]];
        let fi = 0;
        frame.coords.extend(
            (fi * natoms..(fi + 1) * natoms)
                .map(|ai| ai * 3) // Make sure to account for the three values per position.
                .map(|ai| [0.0, 1.0, 2.0].map(|v| ai as f32 + v)),
        );
        frame
    });

    for frame in frames {
        traj.write(&frame).unwrap();
    }

    Ok(())
}
