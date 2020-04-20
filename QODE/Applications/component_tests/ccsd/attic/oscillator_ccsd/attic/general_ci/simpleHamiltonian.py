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
import copy
import numpy as np
from Applications.general_ci import config, state, generate_config
from math import sqrt



def state_dot(state_obj1, state_obj2):
	dot_prod = 0.0
	size = state_obj1.num_config
	# print("vector size =", size)
	for i in range(size):
		dot_prod += state_obj1.coeffs[i] * state_obj2.coeffs[i]
	return dot_prod


def one_elec_hamiltonian(state_obj, h_mat):
	#
	# h1 = \sum_pq h_pq p+ q
	#
	#    h_pq = Kinetic E + Nuclear Attraction E
	#
	num_spin_orb = h_mat.shape[0]
	# print("num spin orb =", num_spin_orb)
	new_state = state.state(state_obj.get_configs())
	#
	#
	for p in range(num_spin_orb):
		for q in range(num_spin_orb):
			op_list    = [ [config.create, p],[config.annihilate, q] ]
			temp_state = state.operate(op_list, state_obj)
			new_state.add_to( h_mat[p,q], temp_state )
	return new_state


def two_elec_hamiltonian(state_obj, V_mat):
	#
	# h2 = \sum_pqrs 1/4 V_pqrs p+ q+ r s
	#
	#    V_pqrs = V_mat[pqsr] - V_mat[pqrs]
	#
	num_spin_orb = int( sqrt(V_mat.shape[0]) )
	# print("num spin orb =", num_spin_orb)
	new_state = state.state(state_obj.get_configs())
	#
	#
	for p in range(num_spin_orb):
		for q in range(num_spin_orb):
			for r in range(num_spin_orb):
				for s in range(num_spin_orb):
					op_list    = [ [config.create, p], [config.create, q], [config.annihilate, r],[config.annihilate, s] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[p * num_spin_orb + q, s * num_spin_orb + r] - V_mat[p * num_spin_orb + q, r * num_spin_orb + s] 
					new_state.add_to( 0.25 * V_value , temp_state )
	return new_state


def hamiltonian(state_obj, h_mat, V_mat):
	new_state = state.state(state_obj.get_configs())
	new_state = state.add_state( new_state, one_elec_hamiltonian(state_obj, h_mat) )
	new_state = state.add_state( new_state, two_elec_hamiltonian(state_obj, V_mat) )
	return new_state



if __name__ == "__main__":
	state_config = generate_config.generate_CISD_configs(1, 1, 8, 'CISD')
	state_obj    = state.state(state_config)
	#
	#
	#  Now Do CISD Code.
	#
	num_config = len(state_config)
	#
	#
	h_mat = np.load('h_mat.npy')
	V_mat = np.load('V_mat.npy')
	#
	#
	H_2nd = []
	#
	#
	# LOOP OVER H |j> first and it's a lot faster.
	#
	for j in range(num_config):
		state_j = copy.deepcopy(state_obj)
		state_j.coeffs[j] = 1.0
		# state_j.print_info()
		#
		new_state_j = hamiltonian(state_j, h_mat, V_mat)
		# print("H|%d> = " %(j) )
		# new_state_j.print_info()
		#
		#
		# Now Loop over <i|
		for i in range(num_config):
			state_i = copy.deepcopy(state_obj)
			state_i.coeffs[i] = 1.0
			#
			#
			H_2nd += [ state_dot( state_i, new_state_j ) ]
			print("< %d | H | %d > = %f" %( i,j, H_2nd[-1] ) )


	print(H_2nd)
	H_2nd = np.array(H_2nd)
	H_2nd = np.reshape(H_2nd, (num_config,num_config))
	H_2nd = np.matrix(H_2nd)
	#
	#
	np.save('SecondQuantizedExplicitState_H_mat.npy', H_2nd)
	#
	#
	E,U = np.linalg.eigh(H_2nd)
	print("Lowest Eigen Value =", E[0])

























