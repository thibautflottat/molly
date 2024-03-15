import time

import molly
import numpy as np

# path = "../../trajectories/md.xtc"
path = "../../trajectories/md.xtc"
reader = molly.XTCReader(path)
arr = np.zeros((1000, 1000, 3), dtype=np.float32)

# start = time.time()
# for _ in range(10000):
#     reader.read_frame()
# duration = time.time() - start
# print(f"reading '{path}' took {duration:.6} s")

frames = reader.read_frames(slice(0, 1000, 3), [i for i in range(100)])
natoms = len(frames[0].positions)
frame = frames[0]
print(frame.box)
print(f"read {len(frames)} frames with {natoms} atoms each from {path}")
