import numpy

neutral  = (numpy.load("data/atomic_states/Z_2e.npy").T)[:2,:]
cationic = (numpy.load("data/atomic_states/Z_1e.npy").T)[:3,:]
anionic  = (numpy.load("data/atomic_states/Z_3e.npy").T)[:2,:]
