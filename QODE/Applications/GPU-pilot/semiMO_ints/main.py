# Usage (within a Psi4 conda environment):  python main.py <distance>

import os
import sys
import numpy
import psi4
import qode.atoms.integrals.spatial_to_spin as spatial_to_spin
import qode.atoms.integrals.external_engines.psi4_ints as integrals
from qode.many_body.self_consistent_field.fermionic import RHF_RoothanHall_Nonorthogonal, RHF_RoothanHall_BiorthogonalS, RHF_RoothanHall_Biorthogonal, UHF_RoothanHall_Biorthogonal

def biorthogonalize(H, V, S):
    S_inv = numpy.linalg.inv(S) 

    I = S_inv @ S
    refine = numpy.linalg.inv(I)
    S_inv = refine @ S_inv

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

# Build molecules and get integrals

#n_elec = 4
#Be = """\
#Be
#"""
#S, T, U, V, X = integrals.AO_ints(Be, "6-31G")
#H = T + U

n_elec_c = 8
Be2 = """\
Be
Be  1  {distance:f}
""".format(distance=float(sys.argv[1]))
S_c, T_c, U_c, V_c, X_c = integrals.AO_ints(Be2, "6-31G")
H_c = T_c + U_c
Enuc_c = X_c.mol.nuclear_repulsion_energy()

# Normal AO SCF 1Be to get atomic coeff matrix
#_, _, Be_C = RHF_RoothanHall_Nonorthogonal(n_elec, (S,H,V), thresh=1e-12)
#numpy.save("coeffs.npy", Be_C)
Be_C = numpy.load("atomic_MO_coeffs.npy")	# The commented procedure does not work b/c degenerate MOs mixed relative to states that come out of atomic_states code

# Psi4 energy of chain for reference
psi4.set_output_file("output.dat")
psi4.set_options({"scf_type": "pk", "PRINT_MOS":"True"})
print("psi4_energy of Be chain    = ", psi4.energy("SCF/6-31G", molecule=X_c.mol))

# Normal AO SCF Be chain
Be_c_energy, _, _ = RHF_RoothanHall_Nonorthogonal(n_elec_c, (S_c,H_c,V_c), thresh=1e-12)
print("Nonorthogonal energy       = ", Be_c_energy + Enuc_c)

# Biorthogonal AO SCF Be chain
H_c_bi, V_c_bi = biorthogonalize(H_c, V_c, S_c)
Be_BIenergy, _, _ = RHF_RoothanHall_BiorthogonalS(n_elec_c, (S_c,H_c_bi,V_c_bi), thresh=1e-12)
print("Biorthogonal energy AO     = ", Be_BIenergy + Enuc_c)

# Biorthogonal semi-MO SCF Be chain
S_semiMO, H_semiMO, V_semiMO = to_semiMO_basis(S_c, H_c, V_c, Be_C)
H_semiMO_bi, V_semiMO_bi = biorthogonalize(H_semiMO, V_semiMO, S_semiMO)
E_total_semiMO, _, _ = RHF_RoothanHall_BiorthogonalS(n_elec_c, (S_semiMO,H_semiMO_bi,V_semiMO_bi), thresh=1e-12)
print("Biorthogonal energy semiMO = ",E_total_semiMO + Enuc_c)

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
### back to the show
E_test, _, _ = UHF_RoothanHall_Biorthogonal(n_elec_c, (H_semiMO_bi,V_semiMO_bi), thresh=1e-12)
print("Biorthogonal energy semiMO =", E_test + Enuc_c, "(true biorthogonal UHF)")

# Dump integrals for X-CC (in spin orbitals, alpha block before beta block)
os.mkdir("{}_A".format(sys.argv[1]))
numpy.save("{}_A/biorthogonal_semiMO_631G_Be2_h.npy".format(sys.argv[1]), H_semiMO_bi)
numpy.save("{}_A/biorthogonal_semiMO_631G_Be2_V.npy".format(sys.argv[1]), V_semiMO_bi)
