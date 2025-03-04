import MDAnalysis as MDA
import time
import molly


def setup_readers(path):
    mda_reader = MDA.coordinates.XTC.XTCReader(
        path, convert_units=False, refresh_offsets=True
    )
    molly_reader = molly.XTCReader(path)

    return mda_reader, molly_reader


def read_all(path):
    """Read and verify all frames, frame-by-frame."""

    mda_reader, molly_reader = setup_readers(path)

    mda_frames = []
    molly_frames = []

    for i in range(mda_reader.n_frames):
        start = time.time()
        mda_positions = mda_reader.trajectory[i].positions
        if i % 10 == 0:
            print(f"{i:>5} mda: {(time.time() - start) * 1000:6.03f} ms/frame", end="")

        start = time.time()
        molly_positions = molly_reader.pop_frame().positions
        if i % 10 == 0:
            print(f"     molly: {(time.time() - start) * 1000:6.03f} ms/frame")

        assert mda_positions.tolist() == molly_positions.tolist()
        mda_frames.append(mda_positions.copy())
        molly_frames.append(molly_positions.copy())

    return mda_frames, molly_frames


def read_frames_slice(path, slice):
    _, molly_reader = setup_readers(path)
    mda_frames, _ = read_all(path)
    mda_frames = mda_frames[slice]
    print(f"{len(mda_frames) = }")
    molly_frames = molly_reader.read_frames(frame_selection=slice, atom_selection=None)
    for i, (mda_positions, molly_frame) in enumerate(zip(mda_frames, molly_frames)):
        molly_positions = molly_frame.positions
        print(i)

        assert (
            mda_positions.tolist() == molly_positions.tolist()
        ), f"{mda_positions = }\n{molly_positions = }"


path = "../../tests/trajectories/trajectory_smol.xtc"
read_frames_slice(path, slice(None, None))
read_frames_slice(path, slice(None, 20))
read_frames_slice(path, slice(25, 50))
read_frames_slice(path, slice(None, None, 2))
read_frames_slice(path, slice(None, 20, 2))
read_frames_slice(path, slice(25, 50, 2))
read_frames_slice(path, slice(None, None, 3))
read_frames_slice(path, slice(None, 20, 3))
read_frames_slice(path, slice(25, 50, 3))
