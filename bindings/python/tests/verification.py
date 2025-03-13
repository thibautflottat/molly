import time

import molly
import MDAnalysis as MDA
import numpy as np

# a = np.array([[1, 2], [3, 4]])
# molly.process_array(a)


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
    # molly_frames = molly_reader.read_frames(frame_selection=slice, atom_selection=None)
    molly_frames = molly_reader.read_frames(frame_selection=slice, atom_selection=None)
    help(molly_reader)

    # |  read_into_array(
    # |      self,
    # |      /,
    # |      coordinate_array,
    # |      boxvec_array,
    # |      time_array=None,
    # |      frame_selection=None,
    # |      atom_selection=None
    # |  )
    # |      Read all frames into the provided `numpy.ndarray`.
    # |
    # |      The `coordinate_array` must have a shape of `(nframes, natoms, 3)`.
    # |
    # |      The `boxvec_array` must have a shape of `(nframes, 3, 3)`.
    # |
    # |      Returns `True` if the reading operation was successful.
    # |
    # |      # Note
    # |
    # |      This function can perform the reads in a buffered manner, depending on the value of the
    # |      `buffered` attribute.

    for i, (mda_positions, molly_frame) in enumerate(zip(mda_frames, molly_frames)):
        molly_positions = molly_frame.positions
        print(i)

        assert (
            mda_positions.tolist() == molly_positions.tolist()
        ), f"{mda_positions = }\n{molly_positions = }"


def read_into_array(path, slice):
    _, molly_reader = setup_readers(path)
    mda_frames, _ = read_all(path)
    mda_frames = mda_frames[slice]
    nframes = len(mda_frames)
    natoms = len(mda_frames[0])
    print(f"{nframes = }, {natoms = }")
    molly_frames = np.zeros((nframes, natoms, 3), dtype=np.float32)
    molly_boxvecs = np.zeros((nframes, 3, 3), dtype=np.float32)
    molly_reader.read_into_array(
        molly_frames, molly_boxvecs, frame_selection=slice, atom_selection=None
    )

    for i, (mda_positions, molly_positions) in enumerate(zip(mda_frames, molly_frames)):
        print(i)

        assert (
            mda_positions.tolist() == molly_positions.tolist()
        ), f"{mda_positions = }\n{molly_positions = }"
    print("tada!!")


path = "../../tests/trajectories/trajectory_smol.xtc"
read_into_array(path, slice(None, None))
read_frames_slice(path, slice(None, None))
read_frames_slice(path, slice(None, 20))
read_frames_slice(path, slice(25, 50))
read_frames_slice(path, slice(None, None, 2))
read_frames_slice(path, slice(None, 20, 2))
read_frames_slice(path, slice(25, 50, 2))
read_frames_slice(path, slice(None, None, 3))
read_frames_slice(path, slice(None, 20, 3))
read_frames_slice(path, slice(25, 50, 3))
