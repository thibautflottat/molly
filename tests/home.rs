mod common;
use common::trajectories;

fn home(path: impl AsRef<std::path::Path>) -> std::io::Result<()> {
    let mut reader = molly::XTCReader::open(&path)?;
    let mut frame = molly::Frame::default();

    // Go through the frames a first time.
    assert!(
        reader.read_frame(&mut frame).is_ok(),
        "should read the first frame"
    );
    let mut n1 = 1; // We already read the first frame.
    while reader.read_frame(&mut frame).is_ok() {
        n1 += 1;
    }
    assert!(
        reader.read_frame(&mut frame).is_err(),
        "idiot check, reader should be done by now"
    );

    // "Move along home!"
    reader.home()?;

    // Go through the frames again.
    assert!(
        reader.read_frame(&mut frame).is_ok(),
        "should read the first frame, after going home again"
    );
    let mut n2 = 1; // We already read the first frame.
    while reader.read_frame(&mut frame).is_ok() {
        n2 += 1;
    }

    assert_eq!(n1, n2, "the number of frames that were read should match");

    Ok(())
}

#[test]
fn home_adk() -> std::io::Result<()> {
    home(trajectories::ADK)
}

#[test]
fn home_aux() -> std::io::Result<()> {
    home(trajectories::AUX)
}

#[test]
fn home_cob() -> std::io::Result<()> {
    home(trajectories::COB)
}

#[test]
fn home_smol() -> std::io::Result<()> {
    home(trajectories::SMOL)
}

#[test]
fn home_ten() -> std::io::Result<()> {
    home(trajectories::TEN)
}

#[test]
fn home_xyz() -> std::io::Result<()> {
    home(trajectories::XYZ)
}

#[test]
fn home_bad() -> std::io::Result<()> {
    home(trajectories::BAD)
}

#[test]
fn home_delinyah() -> std::io::Result<()> {
    home(trajectories::DELINYAH)
}
