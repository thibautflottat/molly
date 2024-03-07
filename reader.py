import MDAnalysis as MDA
import sys
import time

start = time.time()
path = sys.argv[1]
trajectory = MDA.coordinates.XTC.XTCReader(path)

for i, t in enumerate(trajectory):
    pass
end = time.time()
delta = end - start
print(delta, "s")
