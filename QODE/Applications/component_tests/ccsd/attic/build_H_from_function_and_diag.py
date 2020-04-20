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
import numpy as np
from qode.fermion_field       import occ_strings
from qode.fermion_field.state import state, dot
from simpleHamiltonian import hamiltonian

# Energy calculated from Q-Chem 4.3:
# SCF energy in the final basis set = -2.8551604262
# CCSD total energy                 = -2.87178917
#
# PSI4 CISD energy                  = -2.8717891505016



occ_orbs = [0,4]
vrt_orbs = [1,2,3, 5,6,7]
cisd_strings = occ_strings.CISD(occ_orbs,vrt_orbs)
num_config = len(cisd_strings)

h_mat = np.load('data/h_mat.npy')
V_mat = np.load('data/V_mat.npy')
H = hamiltonian(h_mat,V_mat,cisd_strings)

# Loop over H|j> first and it's a lot faster
H_mat = np.zeros( shape=(num_config,num_config) )
for j in range(num_config):
	state_j = state(cisd_strings,j)
	H_state_j = H(state_j)
	# Now Loop over <i|
	for i in range(num_config):
		state_i = state(cisd_strings,i)
		H_mat[i,j] = dot( state_i, H_state_j )
		print("< %d | H | %d > = %f" %( i,j, H_mat[i,j] ) )
H_mat = np.matrix(H_mat)
np.save('data/H_mat_from_2ndQuant.npy', H_mat)

E,U = np.linalg.eigh(H_mat)
print("Lowest Eigenvalue =", E[0])
