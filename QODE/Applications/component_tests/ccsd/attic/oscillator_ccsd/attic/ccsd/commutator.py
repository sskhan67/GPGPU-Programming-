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

def T(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs):
	# Take any state_obj but make a new_state
	#
	# Call T1() and T2()
	#  T|CC> = T1|CC> + T2|CC>
	#
	return state.add_state( T1(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs),\
							T2(state_obj, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs) )






def tau(state_obj, list_of_idx, occ_orbs, virt_orbs):
	#
	# len( list_of_idx ) == 2: It's a SINGLE tau = virt_orbs[idx[0]]+ occ_orbs[idx[1]]
	# len( list_of_idx ) == 4: It's a DOUBLE tau = virt_orbs[idx[0]]+ occ_orbs[idx[1]] virt_orbs[idx[2]]+ occ_orbs[idx[3]]
	#
	new_state = state.state(state_obj.get_configs())
	#
	#
	if len(list_of_idx) == 2:
		print("SINGLE EXCITATION")
		op_list    = [ [config.create, list_of_idx[0]], [config.annihilate, list_of_idx[1]]  ]
		temp_state = state.operate(op_list, state_obj)
		new_state.add_to( 1.0, temp_state )
	elif len(list_of_idx) == 4:
		print("DOUBLE EXCITATION")
		op_list    = [[config.create, list_of_idx[0]], [config.annihilate, list_of_idx[1]], [config.create, list_of_idx[2]], [config.annihilate, list_of_idx[3]]]
		temp_state = state.operate(op_list, state_obj)
		new_state.add_to( 1.0, temp_state )
	elif len(list_of_idx) == 0:
		print("GROUND STATE, DO NOTHING, Tau Becomes Trivial.")
		new_state = copy.deepcopy(state_obj)
	else:
		raise ValueError
	return new_state




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
		print("H is an OP List")
		H_op = copy.deepcopy(H)
	else:
		print("H is a primitive function")
		H_op = operator_list( [[H]], [1.0] )

	if isinstance(T, operator_list):
		print("T is an OP List")
		T_op = copy.deepcopy(T)
	else:
		print("T is a primitive function")
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



def omega1_ops(H,T,tau):
	# \Omega_\mu,\nu^(1) (t^(n)) = < \mu | exp(-T^(n)) [H, \tau_\nu] exp(T^(n)) | HF >
	# exp(-T^(n)) H exp(T^(n)) = H + [H,T] + 1/2 [[H,T],T] + 1/6 [[[H,T],T],T] + 1/24 [[[[H,T],T],T],T]
	# Here substitute H with [ H, \tau_\nu ]
	# Return all the permuated operators.
	#
	if isinstance(H, operator_list):
		print("H is an OP List")
		H_op = copy.deepcopy(H)
	else:
		print("H is a primitive function")
		H_op = operator_list( [[H]], [1.0] )

	if isinstance(tau, operator_list):
		print("tau_nu is an OP List")
		tau_op = copy.deepcopy(tau)
	else:
		print("tau is a primitive function")
		tau_op = operator_list( [[tau]], [1.0] )
	# 1) Make [H, \tau_\nu]
	#
	#
	new_H_op = commute(H_op, tau_op)
	# new_H_op.print_info()
	# 2) Plug it into the omega0_ops() 
	#
	#
	#
	return omega0_ops(new_H_op,T)


def ccsd_energy(t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
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
	omega0_list = omega0_ops(simpleHamiltonian.hamiltonian, T)
	hf_state    = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	#
	return_state = state.state(t_amp_obj.get_configs())
	#

	input_list = []
	for i in range(len(omega0_list)):
		input_list += [[ omega0_list[i], t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat ]]
	

	omega0_pool = Pool( len(omega0_list) )
	list_return_states = omega0_pool.map(omega0_wrapper, input_list)

	for i in range(len(list_return_states)):
		return_state = state.add_state( return_state, list_return_states[i])



	return simpleHamiltonian.state_dot(hf_state, return_state)


# def omega0_mat(t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
# 	#
# 	#
# 	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
# 		generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
# 	#
# 	#
# 	occ_orbs  = alpha_orbs + beta_orbs
# 	virt_orbs = alpha_virt_orbs + beta_virt_orbs
# 	#
# 	dimension = get_num_reference() + get_num_singles(num_occ_orb,num_virt_orb) + get_num_doubles(num_occ_orb,num_virt_orb) 
# 	Omega0    = np.matrix( np.zeros((dimension-1,1)) )
# 	#
# 	#
# 	omega0_list = omega0_ops(simpleHamiltonian.hamiltonian, T)
# 	hf_state    = state.state(t_amp_obj.get_configs())
# 	hf_state.coeffs[0] = 1.0
# 	#
# 	#
# 	return_state = state.state(t_amp_obj.get_configs())
# 	#
# 	for op_list_obj in omega0_list:
# 		# Top level operator list element
# 		print("-----------------------")
# 		operators = op_list_obj.get_operators()
# 		op_coeffs = op_list_obj.get_coeffs()
# 		for i in range(len(op_coeffs)):
# 			# Each term in the long expressions, accumulate them into return_state
# 			# print(op_coeffs[i], operators[i])
# 			if operators[i][-1] == simpleHamiltonian.hamiltonian:
# 				print("H")
# 				temp_state = operators[i][-1](hf_state, h_mat, V_mat)
# 			else:
# 				print("T")
# 				temp_state = operators[i][-1](hf_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs)
# 			#
# 			#
# 			#
# 			for j in  reversed(range(len(operators[i])-1)):
# 				# Each operator in current term, accumulate them into temp_state
# 				# Now act on states.
# 				#
# 				if operators[i][j] == simpleHamiltonian.hamiltonian:
# 					print("H")
# 					temp_state = operators[i][j]( temp_state, h_mat, V_mat) 
# 				else:
# 					print("T")
# 					temp_state = operators[i][j]( temp_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs) 
# 				################# DON'T FORGET THE COEFFICENTS IN FRONT OF OPERATOR TERMS !!!!
# 			return_state.add_to(op_coeffs[i], temp_state)

# 	for i in range(dimension-1):
# 		bra_state             = state.state(t_amp_obj.get_configs())
# 		bra_state.coeffs[i+1]   = 1.0
# 		Omega0[i,0] = simpleHamiltonian.state_dot(bra_state, return_state)
# 	#
# 	#
# 	return Omega0


def omega0_op_term_onto_HF(op_list_obj, t_amp_obj ,num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat):
	return_state = state.state(t_amp_obj.get_configs())
	hf_state  = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	operators = op_list_obj.get_operators()
	op_coeffs = op_list_obj.get_coeffs()
	for i in range(len(op_coeffs)):
		# Each term in the long expressions, accumulate them into return_state
		# print(op_coeffs[i], operators[i])
		if operators[i][-1] == simpleHamiltonian.hamiltonian:
			print("H")
			temp_state = operators[i][-1](hf_state, h_mat, V_mat)
		else:
			print("T")
			temp_state = operators[i][-1](hf_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs)
		#
		#
		#
		for j in  reversed(range(len(operators[i])-1)):
			# Each operator in current term, accumulate them into temp_state
			# Now act on states.
			#
			if operators[i][j] == simpleHamiltonian.hamiltonian:
				print("H")
				temp_state = operators[i][j]( temp_state, h_mat, V_mat) 
			else:
				print("T")
				temp_state = operators[i][j]( temp_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs) 
		return_state.add_to(op_coeffs[i], temp_state)
	return return_state


def omega0_wrapper(this_input_list):
	return omega0_op_term_onto_HF( this_input_list[0], this_input_list[1], this_input_list[2], this_input_list[3], \
								this_input_list[4], this_input_list[5], this_input_list[6], this_input_list[7] )

def omega0_mat_parallel(t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
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
	omega0_list = omega0_ops(simpleHamiltonian.hamiltonian, T)
	hf_state    = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	#
	return_state = state.state(t_amp_obj.get_configs())
	#
	# for op_list_obj in omega0_list:
	# 	# Top level operator list element
	# 	print("-----------------------")

	input_list = []
	for i in range(len(omega0_list)):
		input_list += [[ omega0_list[i], t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat ]]
	
	print("NUM TERMS =", len(omega0_list))
	omega0_pool = Pool( len(omega0_list) )
	list_return_states = omega0_pool.map(omega0_wrapper, input_list)

	for i in range(len(list_return_states)):
		return_state = state.add_state( return_state, list_return_states[i])

	for i in range(dimension-1):
		bra_state               = state.state(t_amp_obj.get_configs())
		bra_state.coeffs[i+1]   = 1.0
		Omega0[i,0] = simpleHamiltonian.state_dot(bra_state, return_state)
	#
	#
	return Omega0


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


# def omega1_mat(t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
# 	#
# 	#
# 	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
# 		generate_config.get_occ_virt_orbs( num_alpha_elec, num_beta_elec, num_spin_orb )
# 	#
# 	#
# 	occ_orbs  = alpha_orbs + beta_orbs
# 	virt_orbs = alpha_virt_orbs + beta_virt_orbs
# 	#
# 	dimension = get_num_reference() + get_num_singles(num_occ_orb,num_virt_orb) + get_num_doubles(num_occ_orb,num_virt_orb) 
# 	Omega1    = np.matrix( np.zeros((dimension-1,dimension-1)) )
# 	#
# 	#
# 	omega1_list = omega1_ops(simpleHamiltonian.hamiltonian,T,tau)
# 	hf_state    = state.state(t_amp_obj.get_configs())
# 	hf_state.coeffs[0] = 1.0
# 	#
# 	#
# 	excitation_idx_list = get_excitation_idx(num_occ_orb, num_virt_orb, occ_orbs, virt_orbs) 
# 	#
# 	#
# 	#
# 	for idx_nu in range(dimension-1):
# 		# First Loop Over 2nd Dimension (largely time savings)
# 		#
# 		return_state = state.state(t_amp_obj.get_configs())
# 		#
# 		for op_list_obj in omega1_list:
# 			print("-----------------------")
# 			op_list_obj.print_info()
# 			operators = op_list_obj.get_operators()
# 			op_coeffs = op_list_obj.get_coeffs()
# 			for i in range(len(op_coeffs)):
# 				# Each term in the long expressions, accumulate them into return_state
# 				# print(op_coeffs[i], operators[i])
# 				if operators[i][-1] == simpleHamiltonian.hamiltonian:
# 					print("H")
# 					temp_state = operators[i][-1](hf_state, h_mat, V_mat)
# 				elif operators[i][-1] == T:
# 					print("T")
# 					temp_state = operators[i][-1](hf_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs)
# 				else: # Must be tau
# 					print("Tau")
# 					temp_state =  operators[i][-1](hf_state, excitation_idx_list[idx_nu] , occ_orbs, virt_orbs)
# 				#
# 				#
# 				#
# 				for j in  reversed(range(len(operators[i])-1)):
# 					# Each operator in current term, accumulate them into temp_state
# 					# Now act on states.
# 					#
# 					if operators[i][j] == simpleHamiltonian.hamiltonian:
# 						print("H")
# 						temp_state = operators[i][j]( temp_state, h_mat, V_mat ) 
# 					elif operators[i][j] == T:
# 						print("T")
# 						temp_state = operators[i][j]( temp_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs ) 
# 					else:
# 						print("Tau")
# 						temp_state =  operators[i][j]( temp_state, excitation_idx_list[idx_nu], occ_orbs, virt_orbs)
# 				#
# 				#
# 				return_state.add_to(op_coeffs[i], temp_state)
# 		#
# 		#
# 		#
# 		for idx_mu in range(dimension-1):
# 			bra_state                   = state.state(t_amp_obj.get_configs())
# 			bra_state.coeffs[i+1]       = 1.0
# 			Omega1[idx_mu, idx_nu]      = simpleHamiltonian.state_dot(bra_state, return_state)

# 	return omega1_mat

#  This is just a copy of the function name to optimize the if conditions:
#  def omega0_op_term_onto_HF(op_list_obj, t_amp_obj ,num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat):
#  def omega1_op_term_onto_HF(omega1_list, t_amp_obj, excitation_idx_list, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat):

def omega1_op_term_onto_HF(omega1_list, t_amp_obj, excitation_idx_list, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat):
	return_state = state.state(t_amp_obj.get_configs())
	#
	#
	hf_state    = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	for op_list_obj in omega1_list:
		print("-----------------------")
		# op_list_obj.print_info()
		operators = op_list_obj.get_operators()
		op_coeffs = op_list_obj.get_coeffs()
		for i in range(len(op_coeffs)):
			# Each term in the long expressions, accumulate them into return_state
			# print(op_coeffs[i], operators[i])
			if operators[i][-1] == simpleHamiltonian.hamiltonian:
				print("H")
				temp_state = operators[i][-1](hf_state, h_mat, V_mat)
			elif operators[i][-1] == T:
				print("T")
				temp_state = operators[i][-1](hf_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs)
			else: # Must be tau
				print("Tau")
				temp_state =  operators[i][-1](hf_state, excitation_idx_list, occ_orbs, virt_orbs)
			#
			#
			#
			for j in  reversed(range(len(operators[i])-1)):
				# Each operator in current term, accumulate them into temp_state
				# Now act on states.
				#
				if operators[i][j] == simpleHamiltonian.hamiltonian:
					print("H")
					temp_state = operators[i][j]( temp_state, h_mat, V_mat ) 
				elif operators[i][j] == T:
					print("T")
					temp_state = operators[i][j]( temp_state, t_amp_obj, num_occ_orb, num_virt_orb, occ_orbs, virt_orbs ) 
				else:
					print("Tau")
					temp_state =  operators[i][j]( temp_state, excitation_idx_list, occ_orbs, virt_orbs)
			#
			#
			return_state.add_to(op_coeffs[i], temp_state)
	return return_state

def omega1_wrapper(this_input_list):
	return omega1_op_term_onto_HF( this_input_list[0], this_input_list[1], this_input_list[2], this_input_list[3], \
								this_input_list[4], this_input_list[5], this_input_list[6], this_input_list[7], this_input_list[8]  )



def omega1_mat_parallel(t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat):
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
	omega1_list = omega1_ops(simpleHamiltonian.hamiltonian,T,tau)
	hf_state    = state.state(t_amp_obj.get_configs())
	hf_state.coeffs[0] = 1.0
	#
	#
	excitation_idx_list = get_excitation_idx(num_occ_orb, num_virt_orb, occ_orbs, virt_orbs) 
	#
	#
	input_list = []
	for i in range(dimension-1):
		input_list += [ [omega1_list, t_amp_obj, excitation_idx_list[i], num_occ_orb, num_virt_orb, occ_orbs, virt_orbs, h_mat, V_mat] ]

	#
	#
	omega1_pool = Pool( dimension-1 )
	list_return_states = omega1_pool.map(omega1_wrapper, input_list)
	#
	#
	for idx_mu in range(dimension-1):
		bra_state                   = state.state(t_amp_obj.get_configs())
		bra_state.coeffs[idx_mu+1]  = 1.0
		for idx_nu in range(dimension-1):
			Omega1[idx_mu, idx_nu]    = simpleHamiltonian.state_dot(bra_state, list_return_states[idx_nu])

	return Omega1



if __name__ == "__main__":
	import time
	import numpy as np
	from Applications.general_ci import generate_config, simpleHamiltonian
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

	t0_start = time.time()
	omega0_vector = omega0_mat_parallel(t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
	t0_end   = time.time()
	np.savetxt("OMEGA0_time.txt", np.array([[t0_end - t0_start]]) )

	print(omega0_vector)

	t1_start = time.time()
	omega1_matrix = omega1_mat_parallel(t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
	t1_end   = time.time()
	np.savetxt("OMEGA1_time.txt", np.array([[t1_end - t1_start]]) )

	omega1_matrix_inv = np.linalg.inv(omega1_matrix)

	dt_amps   = -1.0 * omega1_matrix_inv * omega0_vector
	dt_amps_long = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
	t_amp_obj.update_coeffs( t_amp_obj.get_coeffs() + dt_amps_long.A1 )

	E_CCSD = ccsd_energy( t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)

	print("E_CCSD =", E_CCSD)
	cyc_count += 1



	converged = False
	while not converged:
		print("CCSD CYCLE", cyc_count)
		omega0_vector = omega0_mat_parallel(t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
		omega1_matrix = omega1_mat_parallel(t_amp_obj ,num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
		omega1_matrix_inv = np.linalg.inv(omega1_matrix)
		dt_amps       = -1.0 * omega1_matrix_inv * omega0_vector
		dt_amps_long  = np.concatenate( (np.zeros((1,1)), dt_amps), axis=0)
		t_amp_obj.update_coeffs( t_amp_obj.get_coeffs() + dt_amps_long.A1 )
		E_CCSD_new = ccsd_energy( t_amp_obj, num_alpha_elec, num_beta_elec, num_spin_orb, h_mat, V_mat)
		print("E_CCSD =", E_CCSD_new )
		if abs(E_CCSD - E_CCSD_new) /abs(E_CCSD) < 1e-6:
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
	f.write("CCSD converged at Cycle %d" %(cyc_count) )
	f.write("CCSD Energy = %.16f Hartree" %(E_CCSD_new)  )
	f.close()


