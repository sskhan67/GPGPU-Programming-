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
import sys
import numpy as np
from qode.util import parallel
from qode.fermion_field import occ_strings
from qode.fermion_field.state import state, dot
from qode.fermion_field.nested_operator import multiply as op_mult
from hamiltonian import hamiltonian



resources = parallel.resources(int(sys.argv[1]))

occ_orbs = [0,4]
vrt_orbs = [1,2,3, 5,6,7]

cisd_strings = occ_strings.CISD(occ_orbs,vrt_orbs)		# FCI for 2e-, so contraction of matrix rep is exact within basis
num_config = len(cisd_strings)

h_mat = np.load('data/h_mat.npy')
V_mat = np.load('data/V_mat.npy')
H = hamiltonian(0.0, h_mat, V_mat, cisd_strings, occ_orbs, vrt_orbs)



print("matrix multiply")
H2_mat = np.zeros( shape=(num_config,num_config) )
for j in range(num_config):
	state_j = state(cisd_strings,j)
	H2_state_j = H(H(state_j,resources),resources)
	for i in range(num_config):
		state_i = state(cisd_strings,i)
		H2_mat[i,j] = dot( state_i, H2_state_j )
		print("< %d | H^2 | %d > = %f" %( i,j, H2_mat[i,j] ) )
	sys.stdout.flush()
H2_mat = np.matrix(H2_mat)
E,U = np.linalg.eigh(H2_mat)
print("Lowest Eigenvalue =", E[0])

print("operator left multiply")
H2 = op_mult(H,H, resources=resources)
print("product done")
H2_mat = np.zeros( shape=(num_config,num_config) )
for j in range(num_config):
	state_j = state(cisd_strings,j)
	H2_state_j = H2(state_j,resources)
	for i in range(num_config):
		state_i = state(cisd_strings,i)
		H2_mat[i,j] = dot( state_i, H2_state_j )
		print("< %d | H^2 | %d > = %f" %( i,j, H2_mat[i,j] ) )
	sys.stdout.flush()
H2_mat = np.matrix(H2_mat)
E,U = np.linalg.eigh(H2_mat)
print("Lowest Eigenvalue =", E[0])

print("operator right multiply")
H2 = op_mult(H,H, use_right_multiply=True, resources=resources)
print("product done")
H2_mat = np.zeros( shape=(num_config,num_config) )
for j in range(num_config):
	state_j = state(cisd_strings,j)
	H2_state_j = H2(state_j,resources)
	for i in range(num_config):
		state_i = state(cisd_strings,i)
		H2_mat[i,j] = dot( state_i, H2_state_j )
		print("< %d | H^2 | %d > = %f" %( i,j, H2_mat[i,j] ) )
	sys.stdout.flush()
H2_mat = np.matrix(H2_mat)
E,U = np.linalg.eigh(H2_mat)
print("Lowest Eigenvalue =", E[0])
