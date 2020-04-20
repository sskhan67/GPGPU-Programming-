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
import config
import state
import generate_config
from multiprocessing import Pool
#
#
# This is the 14 pieces of the normal order Hamiltonian
#



def ham1(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h1 = E_0
	#
	return state.scaled(state_obj, E_0)



def ham2(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h2 = \sum_ia F_ia i+ a
	#
	# Make a container for the summation
	new_state = state.state(state_obj.get_configs())
	#
	#
	for i in range(num_occ_orb):
		# Get the orbital number 
		orb_i = occ_orbs[i]
		for a in range(num_virt_orb):
			# Get the orbital number
			orb_a = virt_orbs[a]
			#
			# Get things ready
			op_list    = [ [config.create,orb_i], [config.annihilate,orb_a] ]
			temp_state = state.operate(op_list, state_obj)
			new_state.add_to( F_mat[orb_i, orb_a], temp_state )

	return new_state



def ham3(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h3 = \sum_ij -F_ij  j i+
	#
	# Make a container
	new_state = state.state(state_obj.get_configs())
	#
	#
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for j in range(num_occ_orb):
			orb_j = occ_orbs[j]
			#
			#
			op_list    = [ [config.annihilate, orb_j], [config.create, orb_i] ]
			temp_state = state.operate(op_list, state_obj)
			new_state.add_to( -F_mat[orb_i, orb_j], temp_state )
	#
	#
	return new_state



def ham4(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h4 = \sum_ab F_ab a+ b
	#
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for a in range(num_virt_orb):
		orb_a = virt_orbs[a]
		for b in range(num_virt_orb):
			orb_b = virt_orbs[b]
			#
			#
			op_list = [ [config.create,orb_a],[config.annihilate,orb_b] ]
			temp_state = state.operate(op_list, state_obj)
			new_state.add_to( F_mat[orb_a, orb_b], temp_state )
	#
	#
	return new_state



def ham5(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h5 = \sum_ai F_ai a+ i
	#
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for a in range(num_virt_orb):
		orb_a = virt_orbs[a]
		for i in range(num_occ_orb):
			orb_i = occ_orbs[i]
			#
			#
			op_list = [ [config.create, orb_a], [config.annihilate, orb_i] ]
			temp_state = state.operate(op_list, state_obj)
			new_state.add_to( F_mat[orb_a, orb_i], temp_state )
	#
	#
	return new_state



def ham6(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h6 = \sum_ijab 1/4 V_ijab i+ j+ a b
	#
	#  V_ijab = V_mat[ijba] - V_mat[ijab] 
	#
	num_spin_orb = num_occ_orb + num_virt_orb
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for j in range(num_occ_orb):
			orb_j = occ_orbs[j]
			# HOW TO HANDLE if i == j ???????????????????????????
			#
			for a in range(num_virt_orb):
				orb_a = virt_orbs[a]
				for b in range(num_virt_orb):
					orb_b = virt_orbs[b]
					# HOW TO HANDLE if a == b ???????????????????????????
					#
					op_list = [ [config.create, orb_i], [config.create, orb_j], [config.annihilate, orb_a], [config.annihilate, orb_b] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[orb_i*num_spin_orb + orb_j, orb_b*num_spin_orb + orb_a] - V_mat[orb_i*num_spin_orb + orb_j, orb_a*num_spin_orb + orb_b]
					new_state.add_to( 0.25*V_value, temp_state )
	#
	#
	#
	return new_state



def ham7(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h7 = \sum_ijab 1/4 V_abij a+ b+ i j
	#
	#  V_abij = V_mat[abji] - V_mat[abij]
	#
	num_spin_orb = num_occ_orb + num_virt_orb
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for j in range(num_occ_orb):
			orb_j = occ_orbs[j]
			# 
			#
			for a in range(num_virt_orb):
				orb_a = virt_orbs[a]
				for b in range(num_virt_orb):
					orb_b = virt_orbs[b]
					# 
					#
					op_list = [ [config.create, orb_a], [config.create, orb_b], [config.annihilate, orb_i], [config.annihilate, orb_j] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[orb_a*num_spin_orb + orb_b, orb_j*num_spin_orb + orb_i] - V_mat[orb_a*num_spin_orb + orb_b, orb_i*num_spin_orb + orb_j]
					new_state.add_to( 0.25*V_value, temp_state )
	#
	#
	#
	return new_state



def ham8(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h8 = \sum_ijkl 1/4 V_ijkl k l i+ j+
	#
	#  V_ijkl = V_mat[ijlk] - V_mat[ijkl]
	#
	num_spin_orb = num_occ_orb + num_virt_orb
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for j in range(num_occ_orb):
			orb_j = occ_orbs[j]
			# 
			#
			for k in range(num_occ_orb):
				orb_k = occ_orbs[k]
				for l in range(num_occ_orb):
					orb_l = occ_orbs[l]
					# 
					#
					op_list = [ [config.annihilate, orb_k], [config.annihilate, orb_l], [config.create, orb_i], [config.create, orb_j] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[orb_i*num_spin_orb + orb_j, orb_l*num_spin_orb + orb_k] - V_mat[orb_i*num_spin_orb + orb_j, orb_k*num_spin_orb + orb_l]
					new_state.add_to( 0.25*V_value, temp_state )
	#
	#
	#	
	return new_state
	



def ham9(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h9 = \sum_abcd 1/4 V_abcd a+ b+ c d
	#
	#  V_abcd = V_mat[abdc] - V_mat[abcd]
	#
	#
	num_spin_orb = num_occ_orb + num_virt_orb
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for a in range(num_virt_orb):
		orb_a = virt_orbs[a]
		for b in range(num_virt_orb):
			orb_b = virt_orbs[b]
			# 
			#
			for c in range(num_virt_orb):
				orb_c = virt_orbs[c]
				for d in range(num_virt_orb):
					orb_d = virt_orbs[d]
					# 
					#
					op_list = [ [config.create, orb_a], [config.create, orb_b], [config.annihilate, orb_c], [config.annihilate, orb_d] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[orb_a*num_spin_orb + orb_b, orb_d*num_spin_orb + orb_c] - V_mat[orb_a*num_spin_orb + orb_b, orb_c*num_spin_orb + orb_d]
					new_state.add_to( 0.25*V_value, temp_state )
	#
	#
	#
	return new_state



def ham10(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h10 = \sum_ijab 1/4 V_ajbi a+ i j+ b
	#
	#  V_ajbi = V_mat[ajib] - V_mat[ajbi]
	#
	num_spin_orb = num_occ_orb + num_virt_orb
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for j in range(num_occ_orb):
			orb_j = occ_orbs[j]
			# 
			#
			for a in range(num_virt_orb):
				orb_a = virt_orbs[a]
				for b in range(num_virt_orb):
					orb_b = virt_orbs[b]
					# 
					#
					op_list = [ [config.create, orb_a], [config.annihilate, orb_i], [config.create, orb_j], [config.annihilate, orb_b] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[orb_a*num_spin_orb + orb_j, orb_i*num_spin_orb + orb_b] - V_mat[orb_a*num_spin_orb + orb_j, orb_b*num_spin_orb + orb_i]
					new_state.add_to( 0.25*V_value, temp_state )
	#
	#
	#
	return new_state
	



def ham11(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h11 = \sum_ijak 1/2 V_ijak k j+ i+ a
	#
	#  V_ijak = V_mat[ijka] - V_mat[ijak]
	#
	num_spin_orb = num_occ_orb + num_virt_orb
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for j in range(num_occ_orb):
			orb_j = occ_orbs[j]
			#
			#
			for a in range(num_virt_orb):
				orb_a = virt_orbs[a]
				for k in range(num_occ_orb):
					orb_k = occ_orbs[k]
					#
					#
					op_list = [ [config.annihilate, orb_k], [config.create, orb_j], [config.create, orb_i], [config.annihilate, orb_a] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[orb_i*num_spin_orb + orb_j, orb_k*num_spin_orb + orb_a] - V_mat[orb_i*num_spin_orb + orb_j, orb_a*num_spin_orb + orb_k]
					new_state.add_to( 0.5*V_value, temp_state )		
	#
	#
	#
	return new_state
	



def ham12(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h12 = \sum_abci -1/2 V_ibac b+ c i+ a
	#
	#  V_ibac = V_mat[ibca] - V_mat[ibac]
	#
	num_spin_orb = num_occ_orb + num_virt_orb
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for a in range(num_virt_orb):
		orb_a = virt_orbs[a]
		for b in range(num_virt_orb):
			orb_b = virt_orbs[b]
			#
			#
			for c in range(num_virt_orb):
				orb_c = virt_orbs[c]
				for i in range(num_occ_orb):
					orb_i = occ_orbs[i]
					#
					#
					op_list = [ [config.create, orb_b], [config.annihilate, orb_c], [config.create, orb_i], [config.annihilate, orb_a] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[orb_i*num_spin_orb + orb_b, orb_c*num_spin_orb + orb_a] - V_mat[orb_i*num_spin_orb + orb_b, orb_a*num_spin_orb + orb_c]
					new_state.add_to( -0.5*V_value, temp_state )		
	#
	#
	#
	return new_state
	



def ham13(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h13 = \sum_aijk 1/2 V_jaki a+ i k j+
	#
	#  V_jaki = V_mat[jaik] - V_mat[jaki]
	#
	num_spin_orb = num_occ_orb + num_virt_orb
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for a in range(num_virt_orb):
		orb_a = virt_orbs[a]
		for i in range(num_occ_orb):
			orb_i = occ_orbs[i]
			#
			#
			for j in range(num_occ_orb):
				orb_j = occ_orbs[j]
				for k in range(num_occ_orb):
					orb_k = occ_orbs[k]
					#
					#
					op_list = [ [config.create, orb_a], [config.annihilate, orb_i], [config.annihilate, orb_k], [config.create, orb_j] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[orb_j*num_spin_orb + orb_a, orb_i*num_spin_orb + orb_k] - V_mat[orb_j*num_spin_orb + orb_a, orb_k*num_spin_orb + orb_i]
					new_state.add_to( 0.5*V_value, temp_state )	
	#
	#
	#
	return new_state
	



def ham14(state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# h14 = \sum_abci -1/2 V_baci a+ i b+ c
	#
	#  V_baci = V_mat[baic] - V_mat[baci]
	#
	num_spin_orb = num_occ_orb + num_virt_orb
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	for a in range(num_virt_orb):
		orb_a = virt_orbs[a]
		for b in range(num_virt_orb):
			orb_b = virt_orbs[b]
			#
			#
			for c in range(num_virt_orb):
				orb_c = virt_orbs[c]
				for i in range(num_occ_orb):
					orb_i = occ_orbs[i]
					#
					#
					op_list = [ [config.create, orb_a], [config.annihilate, orb_i], [config.create, orb_b], [config.annihilate, orb_c] ]
					temp_state = state.operate(op_list, state_obj)
					V_value = V_mat[orb_b*num_spin_orb + orb_a, orb_i*num_spin_orb + orb_c] - V_mat[orb_b*num_spin_orb + orb_a, orb_c*num_spin_orb + orb_i]
					new_state.add_to( -0.5*V_value, temp_state )	
	#
	#
	#
	return new_state
	








def ham_wrapper(input_list):
	return input_list[0]( input_list[1], input_list[2], input_list[3], input_list[4], input_list[5], input_list[6], input_list[7], input_list[8] )




def hamiltonian(state_obj, E_0, F_mat, V_mat, num_alpha_elec, num_beta_elec, num_spin_orb):
	#
	# 
	# Make an empty copy to store the new state
	#
	new_state = state.state(state_obj.get_configs())
	#
	# All necessary info,
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
			generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
	occ_orbs  = alpha_orbs + beta_orbs
	virt_orbs = alpha_virt_orbs + beta_virt_orbs
	#
	#
	hamils = [ ham1, ham2, ham3, ham4, ham5, ham6, ham7, ham8, ham9, ham10, ham11, ham12, ham13, ham14 ]
	input_var = [ [hamils[i]] + [ state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs ] for i in range(14) ]
	#
	#
	outputs = []
	#
	for i in range(14):
		outputs += [ ham_wrapper(input_var[i]) ]
	#
	#
	#
	for generated_state_i in outputs:
		new_state = state.add_state( new_state, generated_state_i )
	#
	# Return the result
	return new_state




def hamiltonian_multithread(state_obj, E_0, F_mat, V_mat, num_alpha_elec, num_beta_elec, num_spin_orb):
	#
	# 
	# Make an empty copy to store the new state
	#
	new_state = state.state(state_obj.get_configs())
	#
	# All necessary info,
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
			generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
	occ_orbs  = alpha_orbs + beta_orbs
	virt_orbs = alpha_virt_orbs + beta_virt_orbs
	#
	#
	hamils = [ ham1, ham2, ham3, ham4, ham5, ham6, ham7, ham8, ham9, ham10, ham11, ham12, ham13, ham14 ]
	input_var = [ [hamils[i]] + [ state_obj, E_0, F_mat, V_mat, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs ] for i in range(14) ]
	pool = Pool(14)
	outputs = pool.map(ham_wrapper, input_var)
	#
	#
	for generated_state_i in outputs:
		new_state = state.add_state( new_state, generated_state_i )
	#
	# Return the result
	return new_state


# class normal_order_hamiltonian(object):
# 	def __init__(self):
# 		pass






if __name__ == "__main__":
	np.set_printoptions(precision=3,linewidth=270,threshold=np.nan)
	E_0   = -2.8551604262
	F_mat = np.load("F_mat.npy")
	V_mat = np.load("V_mat.npy")
	print(F_mat)
	state_config = generate_config.generate_CISD_configs(1,1,8, 'CISD')
	state_obj    = state.state(state_config)
	state_obj.coeffs[0] = 1.0
	new_state_obj = hamiltonian(state_obj, E_0, F_mat, V_mat, 1, 1, 8)
	new_state_obj.print_info()
	state_obj.coeffs[0] = 0.0
	state_obj.coeffs[1] = 1.0
	new_state_obj = hamiltonian(state_obj, E_0, F_mat, V_mat, 1, 1, 8)
	new_state_obj.print_info()



