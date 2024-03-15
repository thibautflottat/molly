use molly::{Frame, XTCReader};

fn main() -> std::io::Result<()> {
    let path = std::env::args().skip(1).next().expect("path is required");

    let (tx, rx) = std::sync::mpsc::channel::<Frame>();

    let mut reader = XTCReader::open(path)?;
    let _handle = std::thread::spawn(move || -> std::io::Result<()> {
        let mut frame = Frame::default();
        while reader.read_frame(&mut frame).is_ok() {
            tx.send(frame.clone()).expect("should be able to send");
        }

        Ok(())
    });

    for frame in rx {
        eprintln!("got another frame! -> {:?}", frame.boxvec);
    }

    Ok(())
}
