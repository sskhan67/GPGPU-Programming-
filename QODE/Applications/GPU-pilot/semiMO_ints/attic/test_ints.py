# Usage: python3 test_ints.py <integrals directory>

import sys
import numpy

numpy.set_printoptions(linewidth=500,precision=2)

# Integrals in semi-MO spin orbital basis for a dimer (atom A, then atom B)
h1     = numpy.load("{path}/biorthogonal_semiMO_631G_Be2_h.npy".format(path=sys.argv[1]))		# <^p| h |_q>
v2_raw = numpy.load("{path}/biorthogonal_semiMO_631G_Be2_V.npy".format(path=sys.argv[1]))		# <^p| <^q| V |_r> |_s>

n_orb = h1.shape[0]

# For simplicity antisymmetrize

v2 = numpy.zeros((n_orb, n_orb, n_orb, n_orb))
for p in range(n_orb):
	for q in range(n_orb):
		for r in range(n_orb):
			for s in range(n_orb):
				v2[p][q][r][s] =  v2_raw[p][q][r][s] - v2_raw[p][q][s][r]

def repulsion(rho):
	JK = numpy.zeros((n_orb, n_orb))
	for p in range(n_orb):
		for q in range(n_orb):
			for r in range(n_orb):
				for s in range(n_orb):
					JK[p][q] += rho[s][r] * v2[p][r][q][s]
	return JK


occupied = [0, 1, 0+n_orb//4, 1+n_orb//4 ]
rho = numpy.zeros((n_orb, n_orb))
for i in occupied:  rho[i][i] = 1
JK = repulsion(rho)
E = 0
for i in occupied:  E += h1[i][i] + (1/2.)*JK[i][i]
print("energy of one HF atom in field of other nucleus =",E)

occupied = [0, 1, 0+n_orb//4, 1+n_orb//4, 0+n_orb//2, 1+n_orb//2, 0+n_orb//2+n_orb//4, 1+n_orb//2+n_orb//4]
rho = numpy.zeros((n_orb, n_orb))
for i in occupied:  rho[i][i] = 1
JK = repulsion(rho)
E = 0
for i in occupied:  E += h1[i][i] + (1/2.)*JK[i][i]
print("electronic two overlapping  HF atoms (w/o nuclear repulsion) =",E)
