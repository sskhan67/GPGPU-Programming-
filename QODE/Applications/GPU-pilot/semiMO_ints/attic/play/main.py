import sys
import numpy
import psi4
import qode.atoms.integrals.spatial_to_spin as spatial_to_spin
import psi4_ints as integrals

def biorthogonalize(H, V, S):
    S_inv = numpy.linalg.inv(S) 
    H = S_inv @ H
    for _ in range(2):  V = numpy.tensordot(S_inv, V, axes=([1],[1]))	# cycle through first two tensor axes
    return H, V

def to_semiMO_basis(S, H, V, C1):
    n1 = C1.shape[0]
    N  =  S.shape[0] // n1
    C  = numpy.zeros((N*n1, N*n1))
    for i in range(0, N*n1, n1):  C[i:i+n1, i:i+n1] = C1
    S = C.T @ S @ C
    H = C.T @ H @ C
    for _ in range(4):  V = numpy.tensordot(V, C, axes=([0],[0]))	# cycle through the tensor axes (this assumes everything is real)
    return S, H, V

def coords_Be(n, sep):
    coords = ""
    for i in range(n):
        z = i*sep
        coords += "Be  0.0  0.0  {:f}\n".format(z)
    return coords


n = 100
n_orb = 9 # per atom
sep = 4.5

S = integrals.AO_ints(coords_Be(n, sep), "6-31G").overlap()
S_inv = numpy.linalg.inv(S) 

thresh = 1e-8
sparsity = ""
for i in range(n):
    i0 = i*n_orb
    i1 = i0 + n_orb
    row = ""
    for j in range(n):
        j0 = j*n_orb
        j1 = j0 + n_orb
        if numpy.linalg.norm(S_inv[i0:i1, j0:j1],"fro")>thresh:  row += "x "
        else:                                                    row += "  "
    row += "|\n"
    sparsity += row

print(sparsity)

























"""\

Be_C = numpy.load("atomic_MO_coeffs.npy")

# Biorthogonal semi-MO SCF Be chain
S_semiMO, H_semiMO, V_semiMO = to_semiMO_basis(S_c, H_c, V_c, Be_C)
H_semiMO_bi, V_semiMO_bi = biorthogonalize(H_semiMO, V_semiMO, S_semiMO)


# Double-check with undressed UHF method
H_tmp, V_tmp = spatial_to_spin.one_electron(H_semiMO_bi), spatial_to_spin.two_electron(V_semiMO_bi)

### oops, later codes expect integrals blocked primarily by fragment, then (maybe) spin
n_orb = H_tmp.shape[0]
if n_orb%4 != 0:  raise LogicError
chunk = n_orb // 4
H_semiMO_bi = numpy.zeros((n_orb,n_orb))
for p1,p2 in [(0,0),(chunk,2*chunk),(2*chunk,chunk),(3*chunk,3*chunk)]:
	for p in range(chunk):
		for q1,q2 in [(0,0),(chunk,2*chunk),(2*chunk,chunk),(3*chunk,3*chunk)]:
			for q in range(chunk):
				H_semiMO_bi[p1+p, q1+q] = H_tmp[p2+p, q2+q]
V_semiMO_bi = numpy.zeros((n_orb,n_orb,n_orb,n_orb))
for p1,p2 in [(0,0),(chunk,2*chunk),(2*chunk,chunk),(3*chunk,3*chunk)]:
	for p in range(chunk):
		for q1,q2 in [(0,0),(chunk,2*chunk),(2*chunk,chunk),(3*chunk,3*chunk)]:
			for q in range(chunk):
				for r1,r2 in [(0,0),(chunk,2*chunk),(2*chunk,chunk),(3*chunk,3*chunk)]:
					for r in range(chunk):
						for s1,s2 in [(0,0),(chunk,2*chunk),(2*chunk,chunk),(3*chunk,3*chunk)]:
							for s in range(chunk):
								V_semiMO_bi[p1+p, q1+q, r1+r, s1+s] = V_tmp[p2+p, q2+q, r2+r, s2+s]
"""
