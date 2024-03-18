use std::path::Path;

use xdrfile::Trajectory;

mod trajectories {
    pub const ADK: &str = "tests/trajectories/adk_oplsaa.xtc";
    pub const AUX: &str = "tests/trajectories/aux_edr.xtc";
    pub const COB: &str = "tests/trajectories/cobrotoxin.xtc";
    pub const SMOL: &str = "tests/trajectories/trajectory_smol.xtc";
    pub const TEN: &str = "tests/trajectories/xtc_test_only_10_frame_10_atoms.xtc";
    pub const XYZ: &str = "tests/trajectories/xyz_random_walk.xtc";
    pub const BAD: &str = "tests/trajectories/bad.xtc";
}

mod open {
    use super::*;

    fn open(path: impl AsRef<Path>) -> std::io::Result<()> {
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
}

mod home {
    use super::*;

    fn home(path: impl AsRef<Path>) -> std::io::Result<()> {
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
}

mod compare {
    use super::*;

    fn compare(path: impl AsRef<Path>) -> std::io::Result<()> {
        let mut molly_reader = molly::XTCReader::open(&path)?;
        let mut cf_reader = chemfiles::Trajectory::open_with_format(&path, 'r', "XTC")
            .expect("couldn't open file using chemfiles");
        let mut xdr_reader =
            xdrfile::XTCTrajectory::open_read(&path).expect("couldn't open file using xdrfile");

        let mut molly_frame = molly::Frame::default();
        let mut cf_frame = chemfiles::Frame::new();
        let num_atoms = xdr_reader
            .get_num_atoms()
            .expect("couldn't get number of atoms from xdrfile");
        let mut xdr_frame = xdrfile::Frame::with_len(num_atoms);

        while molly_reader.read_frame(&mut molly_frame).is_ok() {
            cf_reader
                .read(&mut cf_frame)
                .expect("couldn't read chemfiles frame");
            xdr_reader
                .read(&mut xdr_frame)
                .expect("couldn't read xdrfile frame");

            let molly_positions = molly_frame
                .coords()
                .map(|c| c.to_array())
                .collect::<Vec<_>>();
            let cf_positions = cf_frame.positions();
            let xdr_positions = xdr_frame.coords.as_slice();

            assert_eq!(
                molly_positions.len(),
                cf_positions.len(),
                "molly and chemfiles read a different number of atoms"
            );
            assert_eq!(molly_positions.len(), xdr_positions.len());
            assert_eq!(cf_positions.len(), xdr_positions.len()); // For clarity.

            dbg!(&molly_positions[..4]); // LEFT OFF (2024-03-17 14:38): It's fucked for the tests/test.xtc. No clue why, and no clue when this crept in.
            for i in 0..molly_positions.len() {
                let molly_pos = molly_positions[i];
                let cf_pos = cf_positions[i].map(|v| (v * 0.1) as f32); // Convert to nm.
                let xdr_pos = xdr_positions[i];

                assert_eq!(
                    molly_pos, cf_pos,
                    "position {i} for molly and chemfiles does not match"
                );
                assert_eq!(
                    molly_pos, xdr_pos,
                    "position {i} for molly and xdrfile does not match"
                );
                assert_eq!(
                    cf_pos, xdr_pos,
                    "position {i} for chemfiles and xdrfile does not match"
                ); // For clarity.
            }
        }

        // Make sure that after molly_reader is done, cf_reader and xdr_reader are also both done.
        assert!(
            molly_reader.read_frame(&mut molly_frame).is_err(),
            "idiot check, molly reader should be done by now"
        );
        assert!(
            cf_reader.read(&mut cf_frame).is_err(),
            "chemfiles reader should be done by now"
        );
        assert!(
            xdr_reader.read(&mut xdr_frame).is_err(),
            "xdrfile reader should be done by now"
        );

        Ok(())
    }

    #[test]
    fn compare_adk() -> std::io::Result<()> {
        compare(trajectories::ADK)
    }

    #[test]
    fn compare_aux() -> std::io::Result<()> {
        compare(trajectories::AUX)
    }

    #[test]
    fn compare_cob() -> std::io::Result<()> {
        compare(trajectories::COB)
    }

    #[test]
    fn compare_smol() -> std::io::Result<()> {
        compare(trajectories::SMOL)
    }

    #[test]
    fn compare_ten() -> std::io::Result<()> {
        compare(trajectories::TEN)
    }

    #[test]
    fn compare_xyz() -> std::io::Result<()> {
        compare(trajectories::XYZ)
    }

    #[test]
    fn compare_bad() -> std::io::Result<()> {
        compare(trajectories::BAD)
    }
}
