import numpy as np
import pysal
import time
import imp
import maxpls4_lck
import sys

np.random.seed(100)
w = pysal.lat2W(20,20)
z = np.random.random_sample((w.n,2))
p = np.ones((w.n,1), float)
floor = 100
maxpls2 = imp.reload(maxpls2)
start_time = time.time()
if __name__ == '__main__':
    solution = maxpls4_lck.Maxp(w, z, floor, floor_variable=p, initial=100)
print("--- %s seconds in total ---" % (time.time() - start_time))
