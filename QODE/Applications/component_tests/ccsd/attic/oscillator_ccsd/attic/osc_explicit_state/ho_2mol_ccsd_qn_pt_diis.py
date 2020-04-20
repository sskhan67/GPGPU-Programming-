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
# System Imports
import copy
import time
import numpy as np
from multiprocessing import Pool

# Oscillator System Imports
from qode.coupled_oscillators import oscillator_system as osys
from qode.coupled_oscillators import hamiltonian_class
from qode.coupled_oscillators.analytical import analytical_solution

# config/state Imports
import Applications.general_ci.config as config
import Applications.general_ci.state  as state

# Generation of CISD Vector Imports
from Applications.osc_explicit_state.generate_ho_config import generate_main, get_lowest_analytical, get_rec_obj
from Applications.osc_explicit_state.generate_ho_config import diag_rec_h0_mat, diag_rec_v_matrix
from Applications.ccsd.ccsd_qn_diis_pt_simple_H import operator_list, omega0_ops, flatten_operator_list

# H/T Operator Imports
from Applications.osc_explicit_state.osc_hamiltonian_operator import hamiltonian as H
from Applications.osc_explicit_state.osc_t_operator import T





#------------------------------ OMEGA 0 PART --------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------



def omega0_op_term_onto_ref(string_ops, ops_coeff, t_amp_obj, h_mat, V_mat, rec_num_states):
	return_state = state.state(t_amp_obj.get_configs())
	hf_state  = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0

	# 
	# def T(state_obj, t_amp_obj, h_mat, V_mat, rec_num_states):
	temp_state = string_ops[-1](hf_state, t_amp_obj, h_mat, V_mat, rec_num_states)
	#
	#
	for i in  reversed(range(len(string_ops)-1)):
		# Each operator in current term, accumulate them into temp_state
		# Now act on states.
		temp_state = string_ops[i]( temp_state, t_amp_obj, h_mat, V_mat, rec_num_states ) 
	#
	#
	return_state.add_to(ops_coeff, temp_state)
	return return_state


def omega0_wrapper(this_input_list):
	return omega0_op_term_onto_ref( *this_input_list )

def omega0_mat_parallel(computed_return_state, t_amp_obj, H, T, h_mat, V_mat, rec_num_states):
	#
	#
	if computed_return_state == None:
		#
		#
		#
		dimension = len(t_amp_obj.get_coeffs()) 
		Omega0    = np.matrix( np.zeros((dimension-1,1)) )
		#
		#
		omega0_list = omega0_ops(H, T)
		flat_omega0_list = flatten_operator_list(omega0_list)
		del omega0_list
		# flat_omega0_list.print_info()
		#
		#
		return_state = state.state(t_amp_obj.get_configs())
		#
		# for op_list_obj in omega0_list:
		# 	# Top level operator list element
		# 	print("-----------------------")
		#
		operators = flat_omega0_list.get_operators()
		coeffs    = flat_omega0_list.get_coeffs()
		#
		input_list = []
		for i in range(len(coeffs)):
			input_list += [[ operators[i], coeffs[i], t_amp_obj, h_mat, V_mat, rec_num_states ]]
		#
		# # print("NUM TERMS =", len(omega0_list))
		omega0_pool = Pool( len(coeffs) )
		list_return_states = omega0_pool.map(omega0_wrapper, input_list)
		omega0_pool.terminate()
		#
		#
		# Combine all states.
		for i in range(len(list_return_states)):
			return_state = state.add_state( return_state, list_return_states[i])
	else:
		return_state = computed_return_state
		dimension    = len(t_amp_obj.get_coeffs()) 
		Omega0       = np.matrix( np.zeros((dimension-1,1)) )
	#
	#
	# Compute Omega_0 Vector
	for i in range(dimension-1):
		bra_state               = state.state(t_amp_obj.get_configs())
		bra_state.coeffs[i+1]   = 1.0
		Omega0[i,0] = state.state_dot(bra_state, return_state)
	#
	#
	return Omega0



#------------------------------ CCSD ENERGY FUNCTION ------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------


def ccsd_energy(t_amp_obj , H, T,  h_mat, V_mat, rec_num_states):
	#
	#
	print("Solving CCSD Energy...")
	#
	#
	dimension = len(t_amp_obj.get_coeffs()) 
	#
	#
	omega0_list = omega0_ops(H, T)
	flat_omega0_list = flatten_operator_list(omega0_list)
	del omega0_list
	# flat_omega0_list.print_info()
	#
	#
	hf_state    = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	#
	return_state = state.state(t_amp_obj.get_configs())
	#
	# for op_list_obj in omega0_list:
	# 	# Top level operator list element
	# 	print("-----------------------")

	operators = flat_omega0_list.get_operators()
	coeffs    = flat_omega0_list.get_coeffs()

	input_list = []
	for i in range(len(coeffs)):
		input_list += [[ operators[i], coeffs[i], t_amp_obj, h_mat, V_mat, rec_num_states ]]
	
	# # print("NUM TERMS =", len(omega0_list))
	omega0_pool = Pool( len(coeffs) )
	list_return_states = omega0_pool.map(omega0_wrapper, input_list)
	omega0_pool.terminate()
	#
	#
	# Combine all states.
	for i in range(len(list_return_states)):
		return_state = state.add_state( return_state, list_return_states[i] )


	return state.state_dot(hf_state, return_state), return_state



#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------


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



def compute_delta_t_vec(omega0_vector, orb_energy_diff_list):
	size    = omega0_vector.shape[0]
	dt_amps = np.matrix(np.zeros((size,1)))
	for i in range(size):
		dt_amps[i,0] = -1.0/orb_energy_diff_list[i] * omega0_vector[i,0]
	return dt_amps





def diis_solver(diis_subspace, t_amps, dt_amps, t_amps_list, diis_error_vecs):
	#
	#
	if t_amps_list == None and diis_error_vecs == None:
		new_t_amps_list     = [ copy.deepcopy(t_amps) ]
		new_diis_error_vecs = [ copy.deepcopy(dt_amps) ]
		num_subspace        = 1
	else:
		if len(t_amps_list) != len(diis_error_vecs):
			print("Data Corruption: Data Not Aligned!")
			raise AssertionError
		if len(diis_error_vecs) < diis_subspace:
			new_t_amps_list     = copy.deepcopy(t_amps_list)     + [ copy.deepcopy(t_amps)  ]  # List + one element
			new_diis_error_vecs = copy.deepcopy(diis_error_vecs) + [ copy.deepcopy(dt_amps) ]
		else:
			new_t_amps_list     = copy.deepcopy(t_amps_list[1:])     + [ copy.deepcopy(t_amps) ] 
			new_diis_error_vecs = copy.deepcopy(diis_error_vecs[1:]) + [ copy.deepcopy(dt_amps) ]
		num_subspace = len(new_diis_error_vecs)

		
	# print("DIIS ERROR VECS =")
	# print(new_diis_error_vecs)

	# print("t_amps list=")
	# print(new_t_amps_list)

	 
	print("NUM SUBSPACE =", num_subspace)
	error_mat = -1.0 * np.matrix( np.ones(( num_subspace+1, num_subspace+1 )) )
	error_mat[-1,-1] = 0.0    # Last element set to zero
	for i in range(num_subspace):
		for j in range(num_subspace):
			error_mat[i,j] = new_diis_error_vecs[i].T * new_diis_error_vecs[j]

	# print("ERROR MATRIX")
	# print(error_mat)
	print(error_mat)

	error_mat_inv = np.linalg.inv(error_mat)
	# print("INVERTED ERROR MATRIX")
	# print(error_mat_inv)

	right_hand_vec = np.matrix( np.zeros((num_subspace+1,1)) )
	right_hand_vec[-1,0]  = -1.0

	# print("Right Hand Side Vec =")
	# print(right_hand_vec)


	weights = (error_mat_inv * right_hand_vec)[:-1,:]  # DO NOT INCLUDE THE LAST ELEMENT WHICH IS LAMBDA

	print("WEIGHTS")
	print(weights)
	
	# Try to normalize the weights but NOT SURE if it is necessary.
	# total_weight = 0.0
	# for i in range(num_subspace):
	# 	total_weight += weights[i,0]
	# norm_weights = copy.deepcopy(weights)
	# norm_weights = norm_weights / total_weight

	# print("normalized weights")
	# print(norm_weights)
	

	t_amps_short = np.matrix( np.zeros(( dt_amps.shape[0], 1)) )

	for i in range(num_subspace):
		t_amps_short += weights[i,0] * ( new_t_amps_list[i] + new_diis_error_vecs[i] ) 

	# print("SHORT t amplitudes =")
	# print(t_amps_short)

	t_amps_long  = np.concatenate( (np.zeros((1,1)), t_amps_short), axis=0)
	# print("LONG t amplitudes =")
	# print(t_amps_long)

	return t_amps_long, new_t_amps_list, new_diis_error_vecs




#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------




def main(h_mat, V_mat, rec_num_states, num_max_subspace=5, num_itr_start_diis=1):
	if num_max_subspace < 1:
		print("DIIS Subspace Must Be Greater or Equal to One")
		raise ValueError
	if num_itr_start_diis < 1:
		print("DIIS Iteration Must At Least Start at Cycle One")
		raise ValueError
	#
	#
	ccsd_diis_subspace = num_max_subspace
	print('==================================================================================')
	print('CCSD Solver Starts ...')
	print("DIIS Max Subspcae Size =", ccsd_diis_subspace)
	print("DIIS Will Be Invoked at Cycle ", num_itr_start_diis)
	#
	#
	#
	# state_config = generate_config.generate_CISD_configs(num_alpha_elec, num_beta_elec, num_spin_orb, 'CISD')
	ho_configs = generate_main( rec_num_states )
	state_obj  = state.state(ho_configs)
	state_obj.coeffs[0] = 1.0
	#
	#
	orb_energy_diff_list = compute_orb_energy_diff(h_mat, rec_num_states)
	print("=============== ORBITAL ENERGY DIFF LIST ===================================")
	print(orb_energy_diff_list)
	#
	#
	cyc_count = 1
	#
	print('==================================================================================')
	print("CCSD CYCLE", cyc_count)
	#
	t_amp_obj =  state.state(state_obj.get_configs())
	#
	#
	# t_ccsd_start = time.time()
	#
	#
	print("Start Solving Omega_0")
	t0_start = time.time()
	omega0_vector = omega0_mat_parallel(None, t_amp_obj, H, T, h_mat, V_mat, rec_num_states)
	t0_end   = time.time()
	# np.savetxt("OMEGA0_time.txt", np.array([[t0_end - t0_start]]) )
	print("Omega_0 Done.")
	#
	# print(omega0_vector)
	#
	diis_error_vecs = []
	# print("Solving Vector dt")
	# omega1_matrix_inv = np.linalg.inv(omega1_matrix)
	dt_amps = compute_delta_t_vec(omega0_vector, orb_energy_diff_list)
	print("dt amplitudes")
	# print(dt_amps)
	#
	#
	if num_itr_start_diis == 1:
		print("Quasi-Newton PT using DIIS ...")
		t_amps = np.matrix(np.zeros((len(ho_configs)-1,1)))
		t_amps_list = None
		diis_error_vecs = None
		t_amps_long, t_amps_list, diis_error_vecs = diis_solver(ccsd_diis_subspace, t_amps, dt_amps, t_amps_list, diis_error_vecs)
		t_amp_obj.update_coeffs( t_amps_long.A1 )
	else:
		print("Quasi-Newton PT Solver")
		dt_amps_long = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
		t_amp_obj.update_coeffs( t_amp_obj.get_coeffs() + dt_amps_long.A1 )

	# t_amp_obj.print_info()
	print("New t Amplitudes Done.")
	E_CCSD, computed_omega0_state = ccsd_energy( t_amp_obj, H, T, h_mat, V_mat, rec_num_states )
	# 
	print('--------------------------------------------------------------')
	print("E_CCSD =", E_CCSD)
	print('--------------------------------------------------------------')
	#
	# 
	#
	cyc_count += 1
	#
	#
	converged = False
	while not converged:
		print('==================================================================================')
		print("CCSD CYCLE", cyc_count)
		print("Start Solving Omega_0")
		omega0_vector = omega0_mat_parallel(computed_omega0_state, t_amp_obj, H, T, h_mat, V_mat, rec_num_states)
		print("Omega_0 Done.")
		# print(omega0_vector)
		#
		#
		dt_amps = compute_delta_t_vec(omega0_vector, orb_energy_diff_list)
		print("dt amplitudes")
		# print(dt_amps)		
		#
		if cyc_count < num_itr_start_diis:
			print("Quasi-Newton PT Solver")
			dt_amps_long = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
			t_amp_obj.update_coeffs( t_amp_obj.get_coeffs() + dt_amps_long.A1 )
		else:
			print("Quasi-Newton PT using DIIS ...")
			#
			#
			if cyc_count == num_itr_start_diis: 
				# Check if t_amps_list and diis_error_vecs they exist
				t_amps = np.matrix(np.zeros((len(ho_configs)-1,1)))
				t_amps_list     = None
				diis_error_vecs = None
				t_amps_long, t_amps_list, diis_error_vecs = diis_solver(ccsd_diis_subspace, t_amps, dt_amps, t_amps_list, diis_error_vecs)
				t_amp_obj.update_coeffs( t_amps_long.A1 )				
			else:
				t_amps_long, t_amps_list, diis_error_vecs = diis_solver(ccsd_diis_subspace, t_amps_long[1:,:], dt_amps, t_amps_list, diis_error_vecs)
				t_amp_obj.update_coeffs( t_amps_long.A1  )
		#
		#
		print("New t Amplitudes Done.")
		E_CCSD_new, computed_omega0_state  = ccsd_energy( t_amp_obj, H, T, h_mat, V_mat, rec_num_states )
		#
		#
		print('--------------------------------------------------------------')
		print("E_CCSD =", E_CCSD_new )
		print('--------------------------------------------------------------')
		#
		#
		E_diff = abs(E_CCSD - E_CCSD_new) /abs(E_CCSD) 
		print("E_diff =", E_diff)
		if E_diff < 1e-12:
			print("CCSD Converged!")
			converged = True
		else:
			E_CCSD = E_CCSD_new
			cyc_count += 1
	
	print("----------------------------------------------------------")	
	print("CCSD converged at Cycle", cyc_count)
	print("CCSD Energy =", E_CCSD_new, 'Hartree')
	print("----------------------------------------------------------")	


	# t_ccsd_end = time.time()
	# np.savetxt("CCSD_TOTAL_TIME.txt", np.array([[ t_ccsd_end - t_ccsd_start ]]) )
	# np.save("ccsd_vec.npy", t_amp_obj.get_coeffs())
	# f = open("ccsd_energy.txt",'w')
	# f.write("CCSD converged at Cycle %d\n" %(cyc_count) )
	# f.write("CCSD Energy = %.16f Hartree" %(E_CCSD_new)  )
	# f.close()
	return E_CCSD_new	




if __name__ == "__main__":
	#
	#
	rec_num_states = [5,5]
	#
	#
	#
	ho1 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho2 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	#
	#
	ho3 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho4 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	#
	#
	coupling_mat = [[0.0, 0.3],[0.3,0.0]]
	#
	#
	mol1 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho1,ho2] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	#
	#
	mol2 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho3,ho4] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	#
	#
	ext_k_mat = [[0.0,0.1],[0.1,0.0]]
	mol_obj = osys.OscillatorSystem( osys.group_of_oscillators( recursive_list =  [mol1, mol2],
                                                                 coupling     =  ext_k_mat ) )
	#
	#
	print('-------------------------------------')
	print("TOTAL k MATRIX")
	for line in mol_obj.get_k_mat():
		print(line)
	print('-------------------------------------')
	print("ANALYTICAL SOLUTION =",get_lowest_analytical(mol_obj))
	print("------------------------------------")
	#
	#
	rec_obj = get_rec_obj( mol_obj, rec_num_states )
	eig1, eig2 = rec_obj.get_sorted_eigen_energy()
	v_mat_list = rec_obj.get_pair_dipole_mat()
	#
	#
	h_mat = diag_rec_h0_mat([eig1, eig2])
	V_mat = diag_rec_v_matrix(rec_num_states, v_mat_list)
	#
	#
	E_CCSD = main(h_mat, V_mat, rec_num_states, num_max_subspace=5, num_itr_start_diis=1)
	
	#
	# np.savetxt("HO_CCSD" + time.ctime() + ".txt", np.matrix([[E_CCSD]]))
	#
	print("------------------------------------")
	print("ANALYTICAL SOLUTION =",get_lowest_analytical(mol_obj))
	print("------------------------------------")
	#

