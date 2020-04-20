#    (C) Copyright 2018 Anthony D. Dutoi
# 
#    This file is part of Qode.
# 
#    Qode is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    Qode is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with Qode.  If not, see <http://www.gnu.org/licenses/>.
#
import pickle
import numpy
from qode.math import *
from qode.atoms.integrals.basis import contractions, basis_function, ON_transformation, transform2, transform4
from qode.atoms.integrals.integrals_better import kinetic, charge, charges, e_repulsion
from qode.fermion_field.occ_strings import all_occ_strings
from qode.fermion_field.state       import configuration, state, create, annihilate, op_string
import CI_Hamiltonian
import recoupled_Hamiltonian
from excite_strings import all_excite_strings
import state as recoupled_state



n_elec_A = 2
n_elec_B = 2

dist     = 0.74 / 0.529177
distance =    4 / 0.529177
rA1 = (0,0,-distance/2-dist/2)
rA2 = (0,0,-distance/2+dist/2)
rB1 = (0,0,+distance/2-dist/2)
rB2 = (0,0,+distance/2+dist/2)

basis = []
basis += [basis_function(contractions["3-21G"]["H 1s"],rA1)]
basis += [basis_function(contractions["3-21G"]["H 2s"],rA1)]
basis += [basis_function(contractions["3-21G"]["H 1s"],rA2)]
basis += [basis_function(contractions["3-21G"]["H 2s"],rA2)]
basis += [basis_function(contractions["3-21G"]["H 1s"],rB1)]
basis += [basis_function(contractions["3-21G"]["H 2s"],rB1)]
basis += [basis_function(contractions["3-21G"]["H 1s"],rB2)]
basis += [basis_function(contractions["3-21G"]["H 2s"],rB2)]
basis_range_A = (0,4)
basis_range_B = (4,8)
basis_A = basis[slice(*basis_range_A)]
basis_B = basis[slice(*basis_range_B)]
n_orb_A = len(basis_A)
n_orb_B = len(basis_B)
n_orb = n_orb_A + n_orb_B

nuclei = []
nuclei += [charge(1,rA1)]
nuclei += [charge(1,rA2)]
nuclei += [charge(1,rB1)]
nuclei += [charge(1,rB2)]

N_A = charges(nuclei[:2])
N_B = charges(nuclei[2:])
T = kinetic()
V = e_repulsion()



T_mat_A = (basis_A|T|basis_A)
N_mat_AA = (basis_A|N_A|basis_A)
N_mat_AB = (basis_A|N_B|basis_A)
h_mat_A = T_mat_A + N_mat_AA
V_mat_A = ((basis_A,basis_A)|V|(basis_A,basis_A))

U_A = ON_transformation(basis_A)
h_mat_A = transform2(h_mat_A, U_A)
V_mat_A = transform4(V_mat_A, U_A)
N_mat_A = transform2(N_mat_AB, U_A)	# cross attraction

occ_strings_A = all_occ_strings(2*n_orb_A, n_elec_A)
n_config_A = len(occ_strings_A)
H_A = CI_Hamiltonian.H(h_mat_A, V_mat_A, occ_strings_A).terms

try:
	mol_A_file = open("mol_A.pickle","rb")
except:
	H_mat_A = numpy.zeros((n_config_A,n_config_A))
	for op,h in H_A:
		for j,J in enumerate(occ_strings_A):
			ket = configuration(J)
			Hket = op | ket
			if Hket.is_not_null():
				occ_string, phase = Hket.value()
				i = occ_strings_A.index(occ_string)
				H_mat_A[i,j] += phase * h

	evals_A,evecs = numpy.linalg.eigh(H_mat_A)
	indices = evals_A.argsort()
	evals_A = evals_A[indices]
	evecs = evecs[:,indices]

	for i in range(n_config_A):  evals_A[i] += N_A.classical_interaction_energy

	eigenstates_A = []
	for j in range(n_config_A):
		eigenstate = state(occ_strings_A)
		for i in range(n_config_A):
			eigenstate.increment(state(occ_strings_A,i), evecs[i,j])
		eigenstates_A += [eigenstate]

	dens_mats_A = []
	for i in range(n_config_A):
		bra = eigenstates_A[i]
		dens_mats_A_i = []
		for j in range(n_config_A):
			ket = eigenstates_A[j]
			dens_mat_A_ij = numpy.zeros((n_orb_A,n_orb_A))
			for p in range(n_orb_A):
				for q in range(n_orb_A):
					op = op_string( create(2*p),   annihilate(2*q)   )	# alpha
					dens_mat_A_ij[p,q] += bra.dot(op|ket)
					op = op_string( create(2*p+1), annihilate(2*q+1) )	# beta
					dens_mat_A_ij[p,q] += bra.dot(op|ket)
			dens_mats_A_i += [dens_mat_A_ij]
		dens_mats_A += [dens_mats_A_i]

	pickle.dump((evals_A,dens_mats_A),open("mol_A.pickle","wb"))
else:
	evals_A,dens_mats_A = pickle.load(mol_A_file)

print(evals_A[0])





T_mat_B = (basis_B|T|basis_B)
N_mat_BB = (basis_B|N_B|basis_B)
N_mat_BA = (basis_B|N_A|basis_B)
h_mat_B = T_mat_B + N_mat_BB
V_mat_B = ((basis_B,basis_B)|V|(basis_B,basis_B))

U_B = ON_transformation(basis_B)
h_mat_B = transform2(h_mat_B, U_B)
V_mat_B = transform4(V_mat_B, U_B)
N_mat_B = transform2(N_mat_BA, U_B)	# cross attraction

occ_strings_B = all_occ_strings(2*n_orb_B, n_elec_B)
n_config_B = len(occ_strings_B)
H_B = CI_Hamiltonian.H(h_mat_B, V_mat_B, occ_strings_B).terms

try:
	mol_B_file = open("mol_B.pickle","rb")
except:
	H_mat_B = numpy.zeros((n_config_B,n_config_B))
	for op,h in H_B:
		for j,J in enumerate(occ_strings_B):
			ket = configuration(J)
			Hket = op | ket
			if Hket.is_not_null():
				occ_string, phase = Hket.value()
				i = occ_strings_B.index(occ_string)
				H_mat_B[i,j] += phase * h

	evals_B,evecs = numpy.linalg.eigh(H_mat_B)
	indices = evals_B.argsort()
	evals_B = evals_B[indices]
	evecs = evecs[:,indices]

	for i in range(n_config_B):  evals_B[i] += N_B.classical_interaction_energy

	eigenstates_B = []
	for j in range(n_config_B):
		eigenstate = state(occ_strings_B)
		for i in range(n_config_B):
			eigenstate.increment(state(occ_strings_B,i), evecs[i,j])
		eigenstates_B += [eigenstate]

	dens_mats_B = []
	for i in range(n_config_B):
		bra = eigenstates_B[i]
		dens_mats_B_i = []
		for j in range(n_config_B):
			ket = eigenstates_B[j]
			dens_mat_B_ij = numpy.zeros((n_orb_B,n_orb_B))
			for p in range(n_orb_B):
				for q in range(n_orb_B):
					op = op_string( create(2*p),   annihilate(2*q)   )	# alpha
					dens_mat_B_ij[p,q] += bra.dot(op|ket)
					op = op_string( create(2*p+1), annihilate(2*q+1) )	# beta
					dens_mat_B_ij[p,q] += bra.dot(op|ket)
			dens_mats_B_i += [dens_mat_B_ij]
		dens_mats_B += [dens_mats_B_i]

	pickle.dump((evals_B,dens_mats_B),open("mol_B.pickle","wb"))
else:
	evals_B,dens_mats_B = pickle.load(mol_B_file)

print(evals_B[0])



U = numpy.zeros((n_orb,n_orb))
for i,m in enumerate(range(*basis_range_A)):
	for j,n in enumerate(range(*basis_range_A)):
		U[m,n] = U_A[i,j]
for i,m in enumerate(range(*basis_range_B)):
	for j,n in enumerate(range(*basis_range_B)):
		U[m,n] = U_A[i,j]

V_mat = ((basis,basis)|V|(basis,basis))
V_mat = transform4(V_mat, U)



n_states_A = 0
while evals_A[n_states_A]<0:  n_states_A += 1
n_states_B = 0
while evals_B[n_states_B]<0:  n_states_B += 1
n_states_A = 5
n_states_B = 5

evals = [evals_A[:n_states_A], evals_B[:n_states_B]]
transition_rhos = [dens_mats_A, dens_mats_B]
basis_ranges = [basis_range_A, basis_range_B]
N_mats = [N_mat_A, N_mat_B]
excite_strings = all_excite_strings([n_states_A,n_states_B])
n_config = len(excite_strings)
print("config space =", n_config)

H = recoupled_Hamiltonian.H(evals, transition_rhos, basis_ranges, V_mat, N_mats, excite_strings).terms

H_mat = numpy.zeros((n_config,n_config))
for op,h in H:
	for j,J in enumerate(excite_strings):
		ket = recoupled_state.configuration(J)
		Hket = op | ket
		if Hket.is_not_null():
			excite_string, phase = Hket.value()
			i = excite_strings.index(excite_string)
			H_mat[i,j] += phase * h

evals,evecs = numpy.linalg.eigh(H_mat)
indices = evals.argsort()
evals = evals[indices]
evecs = evecs[:,indices]

E = evals[0] + charges(nuclei).classical_interaction_energy - N_A.classical_interaction_energy - N_B.classical_interaction_energy

print(E)
print(E - evals_A[0] - evals_B[0] )
print("-7.263112e-05")
