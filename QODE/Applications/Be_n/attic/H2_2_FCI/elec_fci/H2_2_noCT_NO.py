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
import numpy
from qode.math import *
from qode.atoms.integrals.basis import contractions, basis_function, symmON_transformation, transform2, transform4
from qode.atoms.integrals.integrals_better import kinetic, charge, charges, e_repulsion
from qode.fermion_field.occ_strings import all_occ_strings
from qode.fermion_field.state       import state, configuration
import CI_Hamiltonian

n_elec = 4

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

nuclei = []
nuclei += [charge(1,rA1)]
nuclei += [charge(1,rA2)]
nuclei += [charge(1,rB1)]
nuclei += [charge(1,rB2)]

N = charges(nuclei)
T = kinetic()
V = e_repulsion()

T_mat = (basis|T|basis)
N_mat = (basis|N|basis)
h_mat = T_mat + N_mat
V_mat = ((basis,basis)|V|(basis,basis))

#U = symmON_transformation(basis)
n_orb = len(basis)
basis_range_A = (0,4)
basis_range_B = (4,8)
basis_A = basis[slice(*basis_range_A)]
basis_B = basis[slice(*basis_range_B)]
U_A = symmON_transformation(basis_A)
U_B = symmON_transformation(basis_B)
U = numpy.zeros((n_orb,n_orb))
for i,m in enumerate(range(*basis_range_A)):
	for j,n in enumerate(range(*basis_range_A)):
		U[m,n] = U_A[i,j]
for i,m in enumerate(range(*basis_range_B)):
	for j,n in enumerate(range(*basis_range_B)):
		U[m,n] = U_A[i,j]

h_mat = transform2(h_mat, U)
V_mat = transform4(V_mat, U)

#occ_strings = all_occ_strings(2*len(basis), n_elec)
occ_strings_A = all_occ_strings((2*len(basis))//2, n_elec//2)
occ_strings_B = all_occ_strings((2*len(basis))//2, n_elec//2)
occ_strings = []
for occA in occ_strings_A:
	for occB in occ_strings_B:
		occ_strings += [occA+occB]
H = CI_Hamiltonian.H(h_mat, V_mat, occ_strings).terms
H_mat = numpy.zeros((len(occ_strings),len(occ_strings)))

for op,h in H:
	for j,J in enumerate(occ_strings):
		ket = configuration(J)
		Hket = op | ket
		try:
			occ_string, phase = Hket.value()
			i = occ_strings.index(occ_string)
		except:
			pass
		else:
			H_mat[i,j] += phase * h

evals,evecs = numpy.linalg.eigh(H_mat)
print(min(evals)+ N.classical_interaction_energy)
