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
from Applications.general_ci import state, config, simpleHamiltonian
import time
import numpy as np
from Applications.general_ci import generate_config
import cProfile, pstats, io

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


def save_timing(filename, a_string):
	f = open(filename, 'a')
	f.write(a_string)
	f.close()


#----------------------------------------------------------------------------------------------------------------
#-------------------------------- T1, T2, TAU FUNCTION DEFINITIONS ---------------------------------------------- 

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

# def T(state_obj, trivial_list_of_idx, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, trivial_h_mat, trivial_V_mat):
# 	# Take any state_obj but make a new_state
# 	#
# 	# Call T1() and T2()
# 	#  T|CC> = T1|CC> + T2|CC>
# 	#
# 	return state.add_state( T1(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs),\
# 							T2(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs) )


def T1_operator(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, trivial_h_mat, trivial_V_mat):
	return T1(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs)

def T2_operator(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, trivial_h_mat, trivial_V_mat):
	return T2(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs)



# def H(state_obj, trivial_list_of_idx, trivial_t_amp_obj, trivial_num_occ_orb, trivial_num_virt_orb, trivial_occ_orbs, trivial_virt_orbs, h_mat, V_mat):
# 	return simpleHamiltonian.hamiltonian(state_obj, h_mat, V_mat)


def tau():
	# This is a symbolic function does nothing but holding a name in the commutation procedure.
	# It's replaced by tau_operator class object later in the code.
	pass



class tau_operator(object):
	def __init__(self, excitation_idx_list):
		self.list_of_idx = copy.deepcopy(excitation_idx_list)
	def __call__(self, state_obj, trivial_t_amp_obj, trivial_num_occ_orb, trivial_num_virt_orb, trivial_occ_orbs, trivial_virt_orbs, trivial_h_mat, trivial_V_mat):
		#
		#
		new_state = state.state(state_obj.get_configs())
		#
		if len(self.list_of_idx) == 2:
			op_list    = [ [config.create, self.list_of_idx[0]], [config.annihilate, self.list_of_idx[1]]  ]
			temp_state = state.operate(op_list, state_obj)
			new_state.add_to( 1.0, temp_state )
		elif len(self.list_of_idx) == 4:
			# print("DOUBLE EXCITATION")
			op_list    = [[config.create, self.list_of_idx[0]], [config.annihilate, self.list_of_idx[1]], [config.create, self.list_of_idx[2]], [config.annihilate, self.list_of_idx[3]]]
			temp_state = state.operate(op_list, state_obj)
			new_state.add_to( 1.0, temp_state )
		elif len(self.list_of_idx) == 0:
			# print("GROUND STATE, DO NOTHING, Tau Becomes Trivial.")
			new_state = copy.deepcopy(state_obj)
		else:
			raise ValueError
		#
		#
		return new_state
	def get_nT(self):
		if len(self.list_of_idx) == 2:
			return 1
		elif len(self.list_of_idx) == 4:
			return 2
		else:
			raise ValueError



#---------------------------------------------------------------------------------------------------------------------------
#----------------------------- H1 AND H2 WRAPPER FUNCTIONS -----------------------------------------------------------------



def H1_operator(state_obj, trivial_t_amp_obj, trivial_num_occ_orb, trivial_num_virt_orb, trivial_occ_orbs, trivial_virt_orbs, h_mat, trivial_V_mat):
	# def one_elec_hamiltonian(state_obj, h_mat):
	return simpleHamiltonian.one_elec_hamiltonian(state_obj, h_mat)



def H2_operator(state_obj, trivial_t_amp_obj, trivial_num_occ_orb, trivial_num_virt_orb, trivial_occ_orbs, trivial_virt_orbs, trivial_h_mat, V_mat):
	# def two_elec_hamiltonian(state_obj, V_mat):
	return simpleHamiltonian.two_elec_hamiltonian(state_obj, V_mat)



#----------------------------------------------------------------------------------------------------------------------------
#----------------------------- OPERATOR CLASS AND FUNCTIONS -----------------------------------------------------------------

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




def check_excitation_range(list_of_ops):
	#
	# Check H1 or H2
	#
	if H1_operator in list_of_ops:
		p = 1
	else:
		p = 2
	#
	# Get n_T
	#
	n_T = 0
	for item in list_of_ops:
		if item == T1_operator:
			n_T += 1
		elif item == T2_operator:
			n_T += 2
		elif isinstance(item, tau_operator):
			n_T += item.get_nT()

	#
	# Get k
	#
	k = len(list_of_ops) - 1
	#
	#
	low_limit  = n_T - p
	high_limit = n_T + p - k
	#
	#
	#
	if low_limit <= high_limit:
		if low_limit > 2 or high_limit < 1:
			return False
		else:
			# print(list_of_ops,"Limits = [%d, %d] n_T = %d  p = %d k = %d" %(low_limit,high_limit,n_T,p,k) )
			return True
	else:
		return False



def get_commuted_operators(H1, H2, T1, T2):
	#
	# H1, H2 can be [H1,tau] and [H2,tau], respectively.
	# T1, T2 have to be primitive functions.
	#
	if isinstance(H1, operator_list) and isinstance(H2, operator_list):
		op_h = operator_list( H1.get_operators() + H2.get_operators(), H1.get_coeffs() + H2.get_coeffs() )
	else:
		op_h = operator_list([[H1],[H2]], [1.0,1.0]) 
	#
	# T1 and T2 are always functions:
	op_t = operator_list([[T1],[T2]], [1.0,1.0]) 
	#
	#
	#
	op1  = commute(op_h, op_t)
	# term_count += len(op1.get_coeffs())
	all_op_terms = operator_list( op_h.get_operators() + op1.get_operators(), op_h.get_coeffs() + op1.get_coeffs() )
	# print("Commute Once")
	# op1.print_info()
	for i in range(3):
		# print("Commute %d times" %(i+2))
		op1 = commute(op1, op_t)
		# term_count += len(op1.get_coeffs())
		all_op_terms.update_operators( all_op_terms.get_operators() + op1.get_operators() )
		new_coeffs = copy.deepcopy(op1.get_coeffs())
		if i == 0:
			new_coeffs = scale_py_list(new_coeffs, 0.5)
		elif i == 1:
			new_coeffs = scale_py_list(new_coeffs, 1.0/6.0)
		elif i == 2:
			new_coeffs = scale_py_list(new_coeffs, 1.0/24.0)
		all_op_terms.update_coeffs( all_op_terms.get_coeffs() + new_coeffs )
		# op1.print_info()
	#
	#
	return all_op_terms




def non_zero_omega0_ops(H1,H2,T1,T2):
	# \Omega_\mu^(0) (t^(n))   = < \mu | exp(-T^(n)) H exp(T^(n)) | HF >
	# exp(-T^(n)) H exp(T^(n)) = H + [H,T] + 1/2 [[H,T],T] + 1/6 [[[H,T],T],T] + 1/24 [[[[H,T],T],T],T]
	# Return all the permuated operators.
	#
	#
	# H1, H2, T1, T2 are primitive functions.
	#
	all_ops_obj = get_commuted_operators(H1, H2, T1, T2)
	#
	operators = all_ops_obj.get_operators()
	coeffs    = all_ops_obj.get_coeffs()
	non_zero_ops    = []
	non_zero_coeffs = []
	for i in range(len(operators)):
		if check_excitation_range(operators[i]):
			non_zero_ops    += [copy.deepcopy(operators[i])]
			non_zero_coeffs += [coeffs[i]]
	#
	#
	print("TOTAL NUM TERMS =", len(all_ops_obj.get_operators()))
	print("NON-ZERO  TERMS =", len(non_zero_ops))
	#
	#
	return operator_list(non_zero_ops, non_zero_coeffs)





def omega1_ops(H1,H2,T1,T2,tau):
	# \Omega_\mu,\nu^(1) (t^(n)) = < \mu | exp(-T^(n)) [H, \tau_\nu] exp(T^(n)) | HF >
	# exp(-T^(n)) H exp(T^(n)) = H + [H,T] + 1/2 [[H,T],T] + 1/6 [[[H,T],T],T] + 1/24 [[[[H,T],T],T],T]
	# Here substitute H with [ H, \tau_\nu ]
	# Return all the permuated operators.
	#
	# H1, H2, T1, T2, and tau are ALL primitive functions.
	#
	#
	H1_new = operator_list( [[H1,tau], [tau,H1]], [1.0,-1.0] )
	H2_new = operator_list( [[H2,tau], [tau,H2]], [1.0,-1.0] )
	#
	# 
	all_ops_obj = get_commuted_operators(H1_new, H2_new, T1, T2)
	operators = all_ops_obj.get_operators()
	coeffs    = all_ops_obj.get_coeffs()
	return_ops = []
	return_coeffs = []
	#
	# Remove five T Zero Parts:
	# if length > 5, throw them away.
	for i in range(len(coeffs)):
		if len(operators[i]) < 6:
			return_ops += [ operators[i] ]
			return_coeffs += [coeffs[i]]
	return operator_list(return_ops, return_coeffs)







def omega0_op_term_onto_HF(current_operators, current_coeff, t_amp_obj ,num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat):
	t_start = time.time()
	#
	# This function only deals with one term ( i.e. one sequence ) of operators, no list is involved.
	#
	return_state = state.state(t_amp_obj.get_configs())
	hf_state  = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	# Initialize the Loop Elements.
	temp_state = current_operators[-1](hf_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat)
	#
	# Loop over the rest of the operators.
	for i in reversed(range(len(current_operators)-1)):
		temp_state = current_operators[i]( temp_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat) 
	#
	return_state.add_to(current_coeff, temp_state)
	#
	#
	t_end = time.time()
	op_str = ''
	for item in current_operators:
		op_str += str(item)
	save_timing('omega0_timings.txt', 'time = %f for %s\n' %(t_end - t_start, op_str) )
	return return_state


def omega0_wrapper(this_input_list):
	return omega0_op_term_onto_HF( *this_input_list )

def omega0_mat_parallel(t_amp_obj, H1, H2, T1, T2 ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
	#
	#
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
		generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
	#
	#
	occ_orbs  = alpha_orbs + beta_orbs
	virt_orbs = alpha_virt_orbs + beta_virt_orbs
	#
	dimension = get_num_reference() + get_num_singles(num_occ_orb,num_virt_orb) + get_num_doubles(num_occ_orb,num_virt_orb) 
	Omega0    = np.matrix( np.zeros((dimension-1,1)) )
	#
	#
	omega0_all_ops = non_zero_omega0_ops(H1,H2,T1,T2)
	# omega0_all_ops = get_commuted_operators(H1, H2, T1, T2)
	# omega0_all_ops.print_info()

	operators_non_zero  = omega0_all_ops.get_operators()
	coeffs_non_zero     = omega0_all_ops.get_coeffs()
	#
	#
	hf_state    = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	#
	#
	#
	input_list = []
	for i in range(len(coeffs_non_zero)):
		input_list += [[ operators_non_zero[i], coeffs_non_zero[i], t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat ]]
	
	# print("NUM TERMS =", len(coeffs_non_zero))
	#
	#  Set up parallel parameters:
	#
	if len(coeffs_non_zero) > 30:
		print("Omega_0: Running 30 Threads!")
		omega0_pool = Pool( 30 )
	else:
		print("Omega_0: Running %d Threads!" %(len(coeffs_non_zero)) )
		omega0_pool = Pool(len(coeffs_non_zero))
	#
	# Run Parallel Computations:
	#
	#
	list_return_states = omega0_pool.map(omega0_wrapper, input_list)
	#
	# Kill all the pool processes otherwise it will overflow the 1024 limit
	omega0_pool.terminate()
	#
	#
	# Add all results:
	return_state = state.state(t_amp_obj.get_configs())
	#
	for i in range(len(list_return_states)):
		return_state = state.add_state( return_state, list_return_states[i])
	#
	# Acting Bra onto the returned state.
	#
	del omega0_pool
	for i in range(dimension-1):
		bra_state               = state.state(t_amp_obj.get_configs())
		bra_state.coeffs[i+1]   = 1.0
		Omega0[i,0] = simpleHamiltonian.state_dot(bra_state, return_state)
	#
	#
	return Omega0

#-------------------------------------------------------------------------------------------------------
#----------------------------- CCSD ENERGY FUNCTION ----------------------------------------------------

def ccsd_energy(t_amp_obj, H1, H2, T1, T2, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
	#
	#
	#
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
		generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
	#
	#
	occ_orbs  = alpha_orbs + beta_orbs
	virt_orbs = alpha_virt_orbs + beta_virt_orbs
	#
	dimension = get_num_reference() + get_num_singles(num_occ_orb,num_virt_orb) + get_num_doubles(num_occ_orb,num_virt_orb) 
	Omega0    = np.matrix( np.zeros((dimension-1,1)) )
	#
	#
	omega0_all_ops = non_zero_omega0_ops(H1,H2,T1,T2)
	# omega0_all_ops.print_info()
	operators_non_zero  = omega0_all_ops.get_operators()
	coeffs_non_zero     = omega0_all_ops.get_coeffs()
	#
	#
	hf_state    = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	#
	#
	#
	input_list = []
	for i in range(len(coeffs_non_zero)):
		input_list += [[ operators_non_zero[i], coeffs_non_zero[i], t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat ]]
	
	# print("NUM TERMS =", len(coeffs_non_zero))
	#
	#  Set up parallel parameters:
	#
	if len(coeffs_non_zero) > 30:
		print("Omega_0: Running 30 Threads!")
		omega0_pool = Pool( 30 )
	else:
		print("Omega_0: Running %d Threads!" %(len(coeffs_non_zero)) )
		omega0_pool = Pool(len(coeffs_non_zero))
	#
	# Run Parallel Computations:
	#
	# 
	list_return_states = omega0_pool.map(omega0_wrapper, input_list)
	#
	# Kill all the pool processes otherwise it will overflow the 1024 limit
	omega0_pool.terminate()
	#
	#
	del operators_non_zero
	del coeffs_non_zero
	del input_list
	#
	# Add all results:
	return_state = state.state(t_amp_obj.get_configs())
	#
	for i in range(len(list_return_states)):
		return_state = state.add_state( return_state, list_return_states[i])
	del omega0_pool
	#
	#
	return simpleHamiltonian.state_dot(hf_state, return_state)


#-----------------------------------------------------------------------------------------------------------------------
#---------------------------- OMEGA 1 FUNCTIONS ------------------------------------------------------------------------

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





#  This is just a copy of the function name to optimize the if conditions:
#  def omega0_op_term_onto_HF(op_list_obj, t_amp_obj ,num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat):
#  def omega1_op_term_onto_HF(omega1_list, t_amp_obj, excitation_idx_list, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat):

def omega1_op_term_onto_HF(current_operators, current_coeff, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat):
	t_start = time.time()
	return_state = state.state(t_amp_obj.get_configs())
	#
	#
	hf_state    = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	#
	temp_state = current_operators[-1](hf_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat)
	#
	#
	for i in reversed(range(len(current_operators)-1)):
		# Each operator in current term, accumulate them into temp_state
		# Now act on states.
		# (state_obj, trivial_list_of_idx, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, trivial_h_mat, trivial_V_mat)
		temp_state = current_operators[i](temp_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat) 
	#
	#
	return_state.add_to(current_coeff, temp_state)
	t_end = time.time()
	op_str = ''
	for item in current_operators:
		op_str += str(item)
	save_timing('omega1_timings.txt', 'time = %f for %s\n' %(t_end - t_start, op_str) )

	return return_state


def omega1_wrapper(this_input_list):
	return omega1_op_term_onto_HF( *this_input_list )



def omega1_mat_parallel(t_amp_obj, H1, H2, T1, T2, tau, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
	#
	#
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
		generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
	#
	#
	occ_orbs  = alpha_orbs + beta_orbs
	virt_orbs = alpha_virt_orbs + beta_virt_orbs
	#
	dimension = get_num_reference() + get_num_singles(num_occ_orb,num_virt_orb) + get_num_doubles(num_occ_orb,num_virt_orb) 
	Omega1    = np.matrix( np.zeros((dimension-1,dimension-1)) )
	#
	#
	# Commute all operators but treating tau as a symbolic variable.
	#    Only 4 T terms are kept ( 5 T terms are thrown away in omega1_ops() )    
	#
	omega1_all_ops = omega1_ops(H1,H2,T1,T2,tau)
	#

	# This builds all True/False states for Looping over tau.
	excitation_idx_list = get_excitation_idx(num_occ_orb, num_virt_orb, occ_orbs, virt_orbs)
	#
	# Now Looping Over All Configs and generate different lists of operators for each tau.
	#
	operators = omega1_all_ops.get_operators()
	coeffs    = omega1_all_ops.get_coeffs()

	original_count = len(operators) * (dimension-1)
	reduced_count = 0

	full_omega1_ops    = []
	full_omega1_coeffs = []
	for idx_list in excitation_idx_list:
		ops_same_tau = []
		coeffs_same_tau = []
		for i in range(len(coeffs)):
			if tau in operators[i]:
				temp_op = []
				for j in range(len(operators[i])):
					if operators[i][j] != tau:
						temp_op += [ operators[i][j] ]
					else:
						temp_op += [ tau_operator(idx_list) ]
				if check_excitation_range(temp_op):
					#  Check for Non-Zero Operators, throw away zero terms
					ops_same_tau    += [ temp_op ]
					coeffs_same_tau += [coeffs[i]]
			else:
				if check_excitation_range(operators[i]):
					#  Check for Non-Zero Operators, throw away zero terms
					ops_same_tau    += [ operators[i] ]
					coeffs_same_tau += [coeffs[i]]
		#
		#
		full_omega1_ops    += [ ops_same_tau ]
		full_omega1_coeffs += [ coeffs_same_tau ]
		reduced_count += len(full_omega1_coeffs)
	print("TOTAL NUM TERMS =", original_count)
	print("NON-ZERO  TERMS =", reduced_count)
	#
	#
	hf_state    = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	#
	del operators
	del coeffs
	del omega1_all_ops
	#
	input_list = []
	for i in range(len(full_omega1_coeffs)): # Loop over different excitation ( mu )
		temp_list = []
		for j in range(len(full_omega1_coeffs[i])): # Loop over all the operator terms inside current excitation
			temp_list += [ [full_omega1_ops[i][j], full_omega1_coeffs[i][j] , t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat] ]
		input_list += [ temp_list ]
	#
	#
	all_return_states = []
	#
	omega1_pool = Pool(30)
	
	for i in range(len(full_omega1_coeffs)):
		print("Running On Excitation Num %d" %(i))
		# if len(full_omega1_coeffs[i]) > 30:
		# 	print("Running 30 Threads")
		# 	omega1_pool = Pool( 30 )
		# else:
		# 	print("Running %d Threads" %(len(full_omega1_coeffs[i])))
		# 	omega1_pool = Pool( len(full_omega1_coeffs[i]) )
		#
		# Run Parallel Computation!
		#
		# omega1_pool = Pool(30)
		list_return_states = omega1_pool.map(omega1_wrapper, input_list[i])
		# 
		# Kill all the pool processes otherwise it will overflow the 1024 limit
		# omega1_pool.terminate()
		#
		#
		result_state = state.state(t_amp_obj.get_configs())
		for item in list_return_states:
			result_state = state.add_state(result_state, item)
		del list_return_states
		all_return_states += [ result_state ]
	#
	#
	del omega1_pool
	del input_list
	del full_omega1_coeffs
	del full_omega1_ops
	#
	for idx_mu in range(dimension-1):
		bra_state                   = state.state(t_amp_obj.get_configs())
		bra_state.coeffs[idx_mu+1]  = 1.0
		for idx_nu in range(dimension-1):
			Omega1[idx_mu, idx_nu]    = simpleHamiltonian.state_dot(bra_state, all_return_states[idx_nu])

	return Omega1






def main():
	#
	#
	h_mat = np.load("h_mat.npy")
	V_mat = np.load("V_mat.npy")
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
	print("OCC =", occ_orbs)
	print("VIR =", virt_orbs)

	cyc_count = 1

	print("CCSD CYCLE", cyc_count)
	t_amp_obj =  state.state(state_obj.get_configs())


	t_ccsd_start = time.time()

	print("Start Solving Omega_0")
	t0_start = time.time()
	args = [t_amp_obj, H1_operator, H2_operator, T1_operator, T2_operator, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat ]
	prof = cProfile.Profile()
	prof.enable()
	omega0_vector  = prof.runcall(omega0_mat_parallel, *args)
	prof.disable()
	s = io.StringIO()
	sortby = 'cumulative'
	ps = pstats.Stats(prof, stream=s).sort_stats(sortby)
	ps.print_stats()
	print(s.getvalue())
	# prof.dump_stats('omega0_profile.txt')
	# omega0_vector = omega0_mat_parallel(t_amp_obj, H1_operator, H2_operator, T1_operator, T2_operator, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
	t0_end   = time.time()
	np.savetxt("OMEGA0_time.txt", np.array([[t0_end - t0_start]]) )
	print("Omega_0 Done.")

	print(omega0_vector)

	print("Start Solving Omega_1")
	t1_start = time.time()
	args = [t_amp_obj, H1_operator, H2_operator, T1_operator, T2_operator, tau ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat]
	prof = cProfile.Profile()
	prof.enable()
	omega1_matrix = prof.runcall(omega1_mat_parallel, *args)
	# omega1_matrix = omega1_mat_parallel(t_amp_obj, H1_operator, H2_operator, T1_operator, T2_operator, tau ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
	prof.disable()
	s = io.StringIO()
	sortby = 'cumulative'
	ps = pstats.Stats(prof, stream=s).sort_stats(sortby)
	ps.print_stats()
	print(s.getvalue())
	#
	#
	t1_end   = time.time()
	np.savetxt("OMEGA1_time.txt", np.array([[t1_end - t1_start]]) )
	print("Omega_1 Done.")


	print("Solving Vector dt")
	omega1_matrix_inv = np.linalg.inv(omega1_matrix)

	dt_amps   = -1.0 * omega1_matrix_inv * omega0_vector
	dt_amps_long = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
	t_amp_obj.update_coeffs( t_amp_obj.get_coeffs() + dt_amps_long.A1 )
	print("New t Amplitudes Done.")

	E_CCSD = ccsd_energy( t_amp_obj, H1_operator, H2_operator, T1_operator, T2_operator, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat )

	print("E_CCSD =", E_CCSD)
	cyc_count += 1



	converged = False
	while not converged:
		print("CCSD CYCLE", cyc_count)
		print("Start Solving Omega_0")
		omega0_vector = omega0_mat_parallel(t_amp_obj, H1_operator, H2_operator, T1_operator, T2_operator, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
		print("Omega_0 Done.")
		print("Start Solving Omega_1")
		omega1_matrix = omega1_mat_parallel(t_amp_obj, H1_operator, H2_operator, T1_operator, T2_operator, tau, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
		print("Omega_1 Done.")
		print("Solving Vector dt")
		omega1_matrix_inv = np.linalg.inv(omega1_matrix)
		dt_amps       = -1.0 * omega1_matrix_inv * omega0_vector
		dt_amps_long  = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
		t_amp_obj.update_coeffs( t_amp_obj.get_coeffs() + dt_amps_long.A1 )
		print("New t Amplitudes Done.")
		E_CCSD_new = ccsd_energy( t_amp_obj, H1_operator, H2_operator, T1_operator, T2_operator, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat )
		print("E_CCSD =", E_CCSD_new )
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

if __name__ == '__main__':
	main()


