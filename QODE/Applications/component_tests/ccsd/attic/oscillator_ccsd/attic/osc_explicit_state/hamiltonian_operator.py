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
import Applications.general_ci.config as config
import Applications.general_ci.state  as state


def one_elec_hamiltonian(state_obj, h_mat, rec_num_states):
	#
	# h1 = \sum_pq h_pq p+ q
	#
	#    h_pq = Kinetic E + Nuclear Attraction E
	#
	num_mol = len(rec_num_states)
	# print("num spin orb =", num_spin_orb)
	new_state = state.state(state_obj.get_configs())
	#
	#
	idx_shift = 0
	for mol in range(num_mol):
		for p in range(rec_num_states[mol]):
			for q in range(rec_num_states[mol]):
				op_list    = [ [config.create, p+idx_shift],[config.annihilate, q+idx_shift] ]
				temp_state = state.operate(op_list, state_obj)
				new_state.add_to( h_mat(p+idx_shift, q+idx_shift), temp_state )
		idx_shift += rec_num_states[mol]

	return new_state





def two_elec_hamiltonian(state_obj, V_mat, rec_num_states):
	#
	# h2 = \sum_pqrs 1/4 V_pqrs p+ q+ r s
	#
	#    V_pqrs = V_mat[pqsr] - V_mat[pqrs]
	#
	num_mol = len(rec_num_states)
	#
	# print("num spin orb =", num_spin_orb)
	new_state = state.state(state_obj.get_configs())
	#
	#

	first_idx_shift  = 0
	for i in range(num_mol):               # Mol 1
		# Compute the second index shift
		second_idx_shift = 0
		for num in rec_num_states[:i+1]:
			second_idx_shift += num
		#
		for j in range(i+1, num_mol):      # Mol 2
			#
			#
			#
			#
			for p in range(rec_num_states[i]):      # All States in Mol 1
				for q in range(rec_num_states[j]):  # All States in Mol 2
					#
					#
					for r in range(rec_num_states[i]):   # All States in Mol 1
						for s in range(rec_num_states[j]):  # All States in Mol 2
							#
							op_list    = [ [config.create, p+first_idx_shift], [config.create, q+second_idx_shift], \
										   [config.annihilate, r+first_idx_shift],[config.annihilate, s+second_idx_shift] ]
							temp_state = state.operate(op_list, state_obj)
							V_value = V_mat(i, j, p, q, r, s)
							new_state.add_to( V_value , temp_state )
							# print(this_config, first_idx_shift, second_idx_shift, '==>', first_idx_shift+k+1, second_idx_shift+l+1)
			#
			#
			second_idx_shift += rec_num_states[j]
		first_idx_shift += rec_num_states[i]
	#
	#
	#
	return new_state



def hamiltonian(state_obj, t_amp_obj, h_mat, V_mat, rec_num_states):
	new_state = state.state(state_obj.get_configs())
	new_state = state.add_state( new_state, one_elec_hamiltonian(state_obj, h_mat, rec_num_states) )
	new_state = state.add_state( new_state, two_elec_hamiltonian(state_obj, V_mat, rec_num_states) )
	return new_state




def h_wrapper(this_input_list):
	return hamiltonian(*this_input_list)



	