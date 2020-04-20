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
from qode.atoms.integrals.basis import contractions, basis_function, ON_transformation, transform2, transform4
from qode.atoms.integrals.integrals_better import kinetic, charge, charges, e_repulsion



dist = 1. / 0.529177

basis  = []
basis += [basis_function(contractions["STO-3G"]["H"],(0,0,-dist/2))]
basis += [basis_function(contractions["STO-3G"]["H"],(0,0,+dist/2))]

nuclei = []
nuclei += [charge(1,(0,0,-dist/2))]
nuclei += [charge(1,(0,0,+dist/2))]

N = charges(nuclei)
T = kinetic()
V = e_repulsion()

N_mat = (basis|N|basis)
T_mat = (basis|T|basis)
h_mat = T_mat + N_mat
V_mat = ((basis,basis)|V|(basis,basis))

U = ON_transformation(basis)

h_mat = transform2(h_mat, U)
V_mat = transform4(V_mat, U)

E = 2*h_mat[1,1] + V_mat[1,1,1,1] + N.classical_interaction_energy
print(E)
