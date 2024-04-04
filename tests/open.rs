mod common;
use common::trajectories;

fn open(path: impl AsRef<std::path::Path>) -> std::io::Result<()> {
    let _reader = molly::XTCReader::open(&path)?;
    Ok(())
}

#[test]
fn open_adk() -> std::io::Result<()> {
    open(trajectories::ADK)
}

#[test]
fn open_aux() -> std::io::Result<()> {
    open(trajectories::AUX)
}

#[test]
fn open_cob() -> std::io::Result<()> {
    open(trajectories::COB)
}

#[test]
fn open_smol() -> std::io::Result<()> {
    open(trajectories::SMOL)
}

#[test]
fn open_ten() -> std::io::Result<()> {
    open(trajectories::TEN)
}

#[test]
fn open_xyz() -> std::io::Result<()> {
    open(trajectories::XYZ)
}

#[test]
fn open_bad() -> std::io::Result<()> {
    open(trajectories::BAD)
}

#[test]
fn open_delinyah() -> std::io::Result<()> {
    open(trajectories::DELINYAH)
}
