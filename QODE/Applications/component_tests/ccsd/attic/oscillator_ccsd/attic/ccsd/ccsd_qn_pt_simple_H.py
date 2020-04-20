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
from multiprocessing import Pool
import time
from Applications.general_ci import state, config, simpleHamiltonian



def get_num_reference():
	return 1

def get_num_singles(num_occ_orb, num_virt_orb):
	return num_occ_orb * num_virt_orb

def get_num_doubles(num_occ_orb, num_virt_orb):
	num_configs = 1
	for i in range(num_occ_orb):
		num_configs *= (num_virt_orb - i)
	num_configs = num_configs // 2
	return num_configs


def T1(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# Take any state_obj but make a new_state
	#
	# T1 = \sum_{a i} t_ia a+ i
	# 
	#   t_ia are the coeffs in state_obj but only the SINGLE excitation amplitudes.
	#
	new_state = state.state(state_obj.get_configs())
	t_amps    = t_amp_obj.get_coeffs()
	#
	single_idx = get_num_reference()
	#
	#
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		#
		#
		for a in range(num_virt_orb):
			orb_a = virt_orbs[a]
			#
			#
			op_list     = [[config.create, orb_a], [config.annihilate, orb_i]]
			temp_state  = state.operate(op_list, state_obj)
			new_state.add_to( t_amps[ single_idx ] , temp_state )
			single_idx += 1
	#
	#
	return new_state

def T2(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# Take any state_obj but make a new_state
	#
	# T2 = \sum_{a<b, i<j} t_{ijab} a+ i b+ j 
	#
	#   t_{ijab} are the coeffs in state_obj but only the DOUBLE excitation amplitudes.
	#
	new_state = state.state(state_obj.get_configs())
	t_amps    = t_amp_obj.get_coeffs()
	#
	#
	num_reference = get_num_reference()
	num_singles   = get_num_singles(num_occ_orb, num_virt_orb)
	double_idx = num_reference + num_singles
	#
	#
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for j in range(i+1, num_occ_orb):
			orb_j = occ_orbs[j]
			#
			#
			for a in range(num_virt_orb):
				orb_a = virt_orbs[a]
				for b in range(a+1, num_virt_orb):
					orb_b = virt_orbs[b]
					#
					#
					op_list    = [[config.create, orb_a], [config.annihilate, orb_i], [config.create, orb_b], [config.annihilate, orb_j]]
					temp_state = state.operate(op_list, state_obj)
					new_state.add_to( t_amps[ double_idx ], temp_state )
					double_idx += 1
	#
	#
	return new_state




def T(state_obj, trivial_list_of_idx, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, trivial_h_mat, trivial_V_mat):
	# Take any state_obj but make a new_state
	#
	# Call T1() and T2()
	#  T|CC> = T1|CC> + T2|CC>
	#
	return state.add_state( T1(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs),\
							T2(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs) )




def H(state_obj, trivial_list_of_idx, trivial_t_amp_obj, trivial_num_occ_orb, trivial_num_virt_orb, trivial_occ_orbs, trivial_virt_orbs, h_mat, V_mat):
	return simpleHamiltonian.hamiltonian(state_obj, h_mat, V_mat)





class operator_list(object):
	def __init__(self, list_of_operator_list, list_of_coeffs):
		#
		# list_of_operators = [ op1, op2, op3, ... ... ]
		# list_of_coeffs    = [   1,  -1,  -1, ... ... ]
		#
		self.operators = copy.deepcopy( list_of_operator_list )
		self.coeffs    = copy.deepcopy( list_of_coeffs )
		# Sanity check: len of operators == len of coeffs or raise error
		if len(self.operators) != len(self.coeffs):
			print("num of operator list NOT EQUAL to num coeffs")
			raise AssertionError
		
		#
	# def coeff_change_sign(self, op_list_index):
		# self.coeffs[op_list_index]  = -self.coeffs[op_list_index] 
	def get_operators(self):
		return copy.deepcopy(self.operators)
	def get_coeffs(self):
		return copy.deepcopy( self.coeffs )
	def update_operators(self, new_operators):
		self.operators = copy.deepcopy(new_operators)
	def update_coeffs(self, new_coeffs):
		self.coeffs = copy.deepcopy(new_coeffs)
	def combine_coeffs(self):
		new_coeffs    = []
		new_operators = []
		size = len(self.coeffs)
		for i in range(size):
			if self.operators[i] not in new_operators:
				new_operators += [ copy.deepcopy(self.operators[i]) ]
				new_coeffs    += [ copy.deepcopy(self.coeffs[i]) ]
			else:
				j = 0
				checking = True
				new_size = len(new_coeffs)
				while j < new_size and checking:
					if new_operators[j] == self.operators[i]:
						new_coeffs[j] += self.coeffs[i]
						checking = False
					else:
						j += 1
				if j == new_size:
					print("Found in list but failed to add the coeff, horrible error!")
					raise AssertionError
		self.operators = copy.deepcopy(new_operators)
		self.coeffs    = copy.deepcopy(new_coeffs)

	def print_info(self):
		size = len(self.coeffs)
		for i in range(size):
			print("Coeff =",self.coeffs[i], "OPs =", self.operators[i])





def commute( op_list_obj1, op_list_obj2 ):
	#
	# [A, B] = AB - BA
	#  return AB - BA as one operator_list object with a top level list of two element [AB, -BA].
	#
	# print(" ")
	# print("Commute", op_list_obj1, op_list_obj2)
	# print("------------------------------------------------------------------------------------------------")
	#
	first_ops     =  op_list_obj1.get_operators()
	first_coeffs  =  op_list_obj1.get_coeffs()
	second_ops    =  op_list_obj2.get_operators()
	second_coeffs =  op_list_obj2.get_coeffs()
	#
	size1 = len(first_coeffs)
	size2 = len(second_coeffs)
	#
	new_list   = []
	new_coeffs = []
	#
	# AB part:
	for i in range(size1):
		for j in range(size2):
			# print(i,j)
			new_list   += [ first_ops[i] + second_ops[j] ]
			new_coeffs += [ first_coeffs[i] * second_coeffs[j] ]
	#
	# -BA part:
	for i in range(size2):
		for j in range(size1):
			# print(i,j)
			new_list   += [ second_ops[i] + first_ops[j] ]
			new_coeffs += [ -second_coeffs[i] * first_coeffs[j] ]
	#
	#	
	new_op_list_obj = operator_list( new_list, new_coeffs )
	new_op_list_obj.combine_coeffs()
	return new_op_list_obj





def scale_py_list(py_list, value):
	new_list = [ item * value for item in py_list ]
	return new_list





def omega0_ops(H,T):
	# \Omega_\mu^(0) (t^(n))   = < \mu | exp(-T^(n)) H exp(T^(n)) | HF >
	# exp(-T^(n)) H exp(T^(n)) = H + [H,T] + 1/2 [[H,T],T] + 1/6 [[[H,T],T],T] + 1/24 [[[[H,T],T],T],T]
	# Return all the permuated operators.
	#
	# Check if H, T are primitive functions.
	#
	if isinstance(H, operator_list):
		# print("H is an OP List")
		H_op = copy.deepcopy(H)
	else:
		# print("H is a primitive function")
		H_op = operator_list( [[H]], [1.0] )

	if isinstance(T, operator_list):
		# print("T is an OP List")
		T_op = copy.deepcopy(T)
	else:
		# print("T is a primitive function")
		T_op = operator_list( [[T]], [1.0] )
	# H_op.print_info()
	# T_op.print_info()
	#
	temp_op_list = H_op
	operators = [ copy.deepcopy(temp_op_list) ]

	# print(temp_op_list)
	# print(operators)
	# temp_op_list.print_info()
	for i  in range(4):
		temp_op_list = commute(temp_op_list, T_op)
		operators += [ copy.deepcopy(temp_op_list) ]

	operators[2].update_coeffs( scale_py_list( operators[2].get_coeffs(), 0.5      ) )
	operators[3].update_coeffs( scale_py_list( operators[3].get_coeffs(), 1.0/6.0  ) )	
	operators[4].update_coeffs( scale_py_list( operators[4].get_coeffs(), 1.0/24.0 ) )
	
	return operators





def flatten_operator_list(commutator_list):
	flat_op_list = []
	flat_coeffs  = []
	#
	#
	for operator_list_obj in commutator_list:
		operator = operator_list_obj.get_operators()
		coeffs   = operator_list_obj.get_coeffs()
		for i in range(len(coeffs)):
			flat_op_list += [ operator[i] ]
			flat_coeffs  += [ coeffs[i]   ]
	return operator_list(flat_op_list, flat_coeffs)










#------------------------------ OMEGA 0 PART --------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------



def omega0_op_term_onto_HF(string_ops, ops_coeff, t_amp_obj ,num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat):
	return_state = state.state(t_amp_obj.get_configs())
	hf_state  = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0

	temp_state = string_ops[-1](hf_state, [], t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat)
	#
	#
	for i in  reversed(range(len(string_ops)-1)):
		# Each operator in current term, accumulate them into temp_state
		# Now act on states.
		temp_state = string_ops[i]( temp_state, [], t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat) 
	#
	#
	return_state.add_to(ops_coeff, temp_state)
	return return_state


def omega0_wrapper(this_input_list):
	return omega0_op_term_onto_HF( *this_input_list )

def omega0_mat_parallel(computed_return_state, t_amp_obj, H, T, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
	#
	#
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
			generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
	if computed_return_state == None:
		#
		#
		occ_orbs  = alpha_orbs + beta_orbs
		virt_orbs = alpha_virt_orbs + beta_virt_orbs
		#
		dimension = get_num_reference() + get_num_singles(num_occ_orb,num_virt_orb) + get_num_doubles(num_occ_orb,num_virt_orb) 
		Omega0    = np.matrix( np.zeros((dimension-1,1)) )
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
			input_list += [[ operators[i], coeffs[i], t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat ]]
		
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
		dimension = get_num_reference() + get_num_singles(num_occ_orb,num_virt_orb) + get_num_doubles(num_occ_orb,num_virt_orb) 
		Omega0    = np.matrix( np.zeros((dimension-1,1)) )
	#
	#
	# Compute Omega_0 Vector
	for i in range(dimension-1):
		bra_state               = state.state(t_amp_obj.get_configs())
		bra_state.coeffs[i+1]   = 1.0
		Omega0[i,0] = simpleHamiltonian.state_dot(bra_state, return_state)
	#
	#
	return Omega0



#------------------------------ CCSD ENERGY FUNCTION ------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------


def ccsd_energy(t_amp_obj , H, T, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
	#
	#
	print("Solving CCSD Energy...")
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
		generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
	#
	#
	occ_orbs  = alpha_orbs + beta_orbs
	virt_orbs = alpha_virt_orbs + beta_virt_orbs
	#
	dimension = get_num_reference() + get_num_singles(num_occ_orb,num_virt_orb) + get_num_doubles(num_occ_orb,num_virt_orb) 
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
		input_list += [[ operators[i], coeffs[i], t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat ]]
	
	# # print("NUM TERMS =", len(omega0_list))
	omega0_pool = Pool( len(coeffs) )
	list_return_states = omega0_pool.map(omega0_wrapper, input_list)
	omega0_pool.terminate()
	#
	#
	# Combine all states.
	for i in range(len(list_return_states)):
		return_state = state.add_state( return_state, list_return_states[i] )


	return simpleHamiltonian.state_dot(hf_state, return_state), return_state




#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------


def get_excitation_idx(num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	excitation = []
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for a in range(num_virt_orb):
			orb_a = virt_orbs[a]
			excitation += [[orb_a,orb_i]]
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for j in range(i+1, num_occ_orb):
			orb_j = occ_orbs[j]
			for a in range(num_virt_orb):
				orb_a = virt_orbs[a]
				for b in range(a+1, num_virt_orb):
					orb_b = virt_orbs[b]
					excitation += [[orb_a, orb_i, orb_b, orb_j]]
	return excitation





def compute_orb_energy_diff(fock_orb_energy, excitation_idx):
	energy_diff_list = []
	for item in excitation_idx:
		if len(item) == 2:
			energy_diff_list += [ fock_orb_energy[item[0]] - fock_orb_energy[item[1]] ]
		elif len(item) == 4:
			energy_diff_list += [ fock_orb_energy[item[0]] - fock_orb_energy[item[1]] + fock_orb_energy[item[2]] - fock_orb_energy[item[3]] ]
	return energy_diff_list



def compute_delta_t_vec(omega0_vector, orb_energy_diff_list):
	size    = omega0_vector.shape[0]
	dt_amps = np.matrix(np.zeros((size,1)))
	for i in range(size):
		dt_amps[i] = -1.0/orb_energy_diff_list[i] * omega0_vector[i]
	return dt_amps




if __name__ == "__main__":
	#
	# Total MP2 Energy = -2.86814953377 Hartree (my code)
	# MP2 energy       = -2.86814954    Hartree (Qchem4.4)
	#
	import time
	import numpy as np
	from Applications.general_ci import generate_config, simpleHamiltonian
	import qode.SCF_diag_full.hf as scf_routine
	# 
	#
	# input_file = '/home/yhliu/Qode-clone/Applications//helium-cisd/he.in'
	# E_0, KE, NAE, NRE, ERE, F_mat, T_mat, N_mat, V_mat, vHF, num_alpha_elec, num_beta_elec = scf_routine.main(input_file, 1e-10)
	# print(F_mat)
	h_mat = np.load("h_mat.npy")
	V_mat = np.load("V_mat.npy")
	E_hf_orb = np.load("F.npy")
	#
	state_config = generate_config.generate_CISD_configs(1,1,8, 'CISD')
	state_obj    = state.state(state_config)
	state_obj.coeffs[0] = 1.0
	#
	#
	num_alpha_elec = 1
	num_beta_elec  = 1
	num_spin_orb   = 8
	#
	#
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
			generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
	occ_orbs  = alpha_orbs + beta_orbs
	virt_orbs = alpha_virt_orbs + beta_virt_orbs
	print("Occupied =", occ_orbs)
	print("Virtual  =", virt_orbs)

	excitation_idx = get_excitation_idx(num_occ_orb, num_virt_orb, occ_orbs, virt_orbs)
	orb_energy_diff_list = compute_orb_energy_diff(E_hf_orb, excitation_idx)


	cyc_count = 1

	print('--------------------------------------------------------------')
	print("CCSD CYCLE", cyc_count)
	t_amp_obj =  state.state(state_obj.get_configs())


	t_ccsd_start = time.time()

	print("Start Solving Omega_0")
	t0_start = time.time()
	omega0_vector = omega0_mat_parallel(None, t_amp_obj, H, T, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
	t0_end   = time.time()
	np.savetxt("OMEGA0_time.txt", np.array([[t0_end - t0_start]]) )
	print("Omega_0 Done.")

	print(omega0_vector)


	# print("Solving Vector dt")
	# omega1_matrix_inv = np.linalg.inv(omega1_matrix)
	dt_amps = compute_delta_t_vec(omega0_vector, orb_energy_diff_list)
	dt_amps_long = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
	t_amp_obj.update_coeffs( t_amp_obj.get_coeffs() + dt_amps_long.A1 )
	print("New t Amplitudes Done.")
	E_CCSD, computed_omega0_state = ccsd_energy( t_amp_obj , H, T, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
	print('--------------------------------------------------------------')
	print("E_CCSD =", E_CCSD)
	print('--------------------------------------------------------------')

	cyc_count += 1

	
	converged = False
	while not converged:
		print('--------------------------------------------------------------')
		print("CCSD CYCLE", cyc_count)
		print("Start Solving Omega_0")
		omega0_vector = omega0_mat_parallel(computed_omega0_state, t_amp_obj, H, T, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
		print("Omega_0 Done.")
		#
		#
		dt_amps       = compute_delta_t_vec(omega0_vector, orb_energy_diff_list)
		dt_amps_long  = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
		t_amp_obj.update_coeffs( t_amp_obj.get_coeffs() + dt_amps_long.A1 )
		print("New t Amplitudes Done.")
		E_CCSD_new, computed_omega0_state  = ccsd_energy( t_amp_obj, H, T, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
		print('--------------------------------------------------------------')
		print("E_CCSD =", E_CCSD_new )
		print('--------------------------------------------------------------')
		#
		#
		if abs(E_CCSD - E_CCSD_new) /abs(E_CCSD) < 1e-6:
			print("CCSD Converged!")
			converged = True
		else:
			E_CCSD = E_CCSD_new
			cyc_count += 1
	
	print("----------------------------------------------------------")	
	print("CCSD converged at Cycle", cyc_count)
	print("CCSD Energy =", E_CCSD_new, 'Hartree')
	print("----------------------------------------------------------")	


	t_ccsd_end = time.time()

	np.savetxt("CCSD_TOTAL_TIME.txt", np.array([[ t_ccsd_end - t_ccsd_start ]]) )

	# np.save("ccsd_vec.npy", t_amp_obj.get_coeffs())
	f = open("ccsd_energy.txt",'w')
	f.write("CCSD converged at Cycle %d\n" %(cyc_count) )
	f.write("CCSD Energy = %.16f Hartree" %(E_CCSD_new)  )
	f.close()
	

