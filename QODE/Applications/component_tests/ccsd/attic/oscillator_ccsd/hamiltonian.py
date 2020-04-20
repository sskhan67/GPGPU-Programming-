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
from qode.fermion_field import state, config




def one_elec_hamiltonian(state_obj, h_mat, rec_num_states):
	#
	# h1 = \sum_pq h_pq p+ q
	#
	#    h_pq = Kinetic E + Nuclear Attraction E
	#
	num_mol = len(rec_num_states)
	# print("num spin orb =", num_spin_orb)
	occ_string = state_obj.value()[0]
	new_state  = state.state(occ_string)
	new_state.update_coeffs( np.zeros((len(occ_string))))
	#
	#
	idx_shift = 0
	for mol in range(num_mol):
		for p in range(rec_num_states[mol]):
			for q in range(rec_num_states[mol]):
				# op_list    = [ [config.create, p+idx_shift],[config.annihilate, q+idx_shift] ]
				# temp_state = state.operate(op_list, state_obj)
				op = state.op_string( state.create(p+idx_shift), state.annihilate(q+idx_shift) )
				new_state.increment( op(state_obj), h_mat(p+idx_shift, q+idx_shift) )
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
	occ_string = state_obj.value()[0]
	new_state  = state.state(occ_string)
	new_state.update_coeffs( np.zeros((len(occ_string))))
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
							# op_list    = [ [config.create, p+first_idx_shift], [config.create, q+second_idx_shift], \
							# 			   [config.annihilate, r+first_idx_shift],[config.annihilate, s+second_idx_shift] ]
							# temp_state = state.operate(op_list, state_obj)
							op = state.op_string( state.create(p+first_idx_shift), state.create(q+second_idx_shift), state.annihilate(r+first_idx_shift) , state.annihilate(s+second_idx_shift))
							V_value = V_mat(i, j, p, q, r, s)
							new_state.increment( op(state_obj), V_value)
							# print(this_config, first_idx_shift, second_idx_shift, '==>', first_idx_shift+k+1, second_idx_shift+l+1)
			#
			#
			second_idx_shift += rec_num_states[j]
		first_idx_shift += rec_num_states[i]
	#
	#
	#
	return new_state



def hamiltonian_func(state_obj, h_mat, V_mat, rec_num_states):
	new_state = one_elec_hamiltonian(state_obj, h_mat, rec_num_states)
	new_state.increment( two_elec_hamiltonian(state_obj, V_mat, rec_num_states) )
	return new_state

def h_wrapper(this_input_list):
	return hamiltonian_func(*this_input_list)


class hamiltonian(object):
	"""hamiltonian class as a wrapper for hamiltonian function"""
	def __init__(self, h_mat, V_mat, rec_num_states):
		self.h_mat = copy.deepcopy(h_mat)
		self.V_mat = copy.deepcopy(V_mat)
		self.rec_num_states = copy.deepcopy(rec_num_states)
	def __call__(self, state_obj):
		return hamiltonian_func(state_obj, self.h_mat, self.V_mat, self.rec_num_states)
	




def compute_orb_energy_diff(h_mat, rec_num_states):
	energy_diff_list = []
	num_mol = len(rec_num_states)
	#
	# Adding Singles: e_ai = e_a - e_i
	idx_shift = 0
	for i in range(num_mol):
		for j in range(1, rec_num_states[i]):
			energy_diff_list += [ h_mat(idx_shift+j, idx_shift+j) - h_mat(idx_shift, idx_shift) ]
		idx_shift += rec_num_states[i]
	#
	# Adding Doubles: e_aibj = e_a - e_i + e_b - e_j
	first_idx_shift = 0
	for i in range(num_mol):
		#
		#
		second_idx_shift = 0
		for num in rec_num_states[:i+1]:
			second_idx_shift += num
		#
		#
		for j in range(i+1, num_mol):
			#
			#
			for a in range(1, rec_num_states[i]):
				for b in range(1, rec_num_states[j]):
					energy_diff_list += [ h_mat(first_idx_shift+a, first_idx_shift+a) - h_mat(first_idx_shift,first_idx_shift) +\
					                      h_mat(second_idx_shift+b,second_idx_shift+b) - h_mat(second_idx_shift,second_idx_shift) ]
			second_idx_shift += rec_num_states[j]
		first_idx_shift += rec_num_states[i]
	#
	#
	return energy_diff_list


def compute_delta_t_vec(omega0_state, orb_energy_diff_list):
	occ_string, coeffs = omega0_state.value()
	# print(omega0_vector)
	# print(coeffs)
	# print("orb E diff size =", len(orb_energy_diff_list))
	omega0_vector = coeffs
	size    = coeffs.shape[0]
	# print("SIZE =", size)
	dt_amps = np.zeros((size))
	for i in range(size-1):
		dt_amps[i+1] = -1.0/orb_energy_diff_list[i] * omega0_vector[i+1]
	# print(dt_amps)
	return_state = state.state(occ_string)
	return_state.update_coeffs(dt_amps)
	# print(return_state.value()[1])
	return return_state


class inv_fock(object):
	"""inv_fock class to aid in the precesures of """
	def __init__(self, h_mat, rec_num_states):
		self.h_mat = copy.deepcopy(h_mat)
		self.rec_num_states = copy.deepcopy(rec_num_states)
		self.orb_energy_diff_list = compute_orb_energy_diff(self.h_mat, self.rec_num_states)
	def __call__(self, omega0_state):
		return compute_delta_t_vec(omega0_state, self.orb_energy_diff_list)




