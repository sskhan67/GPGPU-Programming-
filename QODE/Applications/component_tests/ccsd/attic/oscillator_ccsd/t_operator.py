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
from Applications.component_tests.oscillator_ccsd.old_qode.generate_ho_config import generate_singles



def T1(state_obj, t_amp_obj, rec_num_states):
	# Take any state_obj but make a new_state
	#
	# T1 = \sum_{a i} t_ia a+ i
	# 
	#   t_ia are the coeffs in state_obj but only the SINGLE excitation amplitudes.
	#
	occ_string = state_obj.value()[0]
	new_state  = state.state(occ_string)
	new_state.update_coeffs( np.zeros(( len(occ_string) )) )
	t_amps    = t_amp_obj.value()[1]
	idx_shift = 0
	num_mol   = len(rec_num_states)
	t_idx     = 1
	for i in range(num_mol):
		for j in range(1, rec_num_states[i]):
			# print('mol =',i,'ref =>',j)
			#
			# new way:
			# op = state.op_string( state.create(orb_a), state.annihilate(orb_i) )
			# new_state.increment( op(state_obj), t_amps[single_idx] )
			#
			#
			# op_list     = [[config.create, idx_shift + j ], [config.annihilate, idx_shift]]
			# temp_state  = state.operate(op_list, state_obj)
			op = state.op_string( state.create(idx_shift + j), state.annihilate(idx_shift) )
			new_state.increment( op(state_obj), t_amps[ t_idx ] )
			t_idx += 1
		#		
		idx_shift += rec_num_states[i]

	return new_state




def T2(state_obj, t_amp_obj, rec_num_states):
	# Take any state_obj but make a new_state
	#
	# T2 = \sum_{a<b, i<j} t_{ijab} a+ i b+ j 
	#
	#   t_{ijab} are the coeffs in state_obj but only the DOUBLE excitation amplitudes.
	#
	occ_string = state_obj.value()[0]
	new_state = state.state(occ_string)
	new_state.update_coeffs( np.zeros(( len(occ_string) )) )
	#
	t_amps    = t_amp_obj.value()[1]
	#
	num_mol   = len(rec_num_states)
	#
	#
	t_idx = 1 + len(generate_singles(rec_num_states))
	first_idx_shift = 0
	for i in range(num_mol):   # Mol 1
		#
		#
		second_idx_shift = 0
		for num in rec_num_states[:i+1]:
			second_idx_shift += num
		#
		#
		for j in range(i+1, num_mol): # Mol 2
			#
			#
			for a in range(1, rec_num_states[i]):
				#
				#
				for b in range(1, rec_num_states[j]):
					#
					#
					# op_list    = [[config.create, first_idx_shift + a], [config.annihilate, first_idx_shift], [config.create, second_idx_shift + b], [config.annihilate, second_idx_shift] ]
					# temp_state = state.operate(op_list, state_obj)
					op = state.op_string( state.create(first_idx_shift + a), state.annihilate(first_idx_shift), state.create(second_idx_shift + b) , state.annihilate(second_idx_shift))
					new_state.increment( op(state_obj), t_amps[ t_idx ] )
					t_idx += 1
			#
			#
			second_idx_shift += rec_num_states[j]
		first_idx_shift += rec_num_states[i]
	#
	#
	return new_state




def T(state_obj, t_amp_obj, rec_num_states):
	new_state = T1(state_obj, t_amp_obj, rec_num_states)
	new_state.increment( T2(state_obj, t_amp_obj, rec_num_states) )
	return new_state



class T_operator(object):
	"""A wrapper for function T"""
	def __init__(self, t_amp_obj, rec_num_states):
		self.t_amp_obj = copy.deepcopy(t_amp_obj)
		self.rec_num_states = copy.deepcopy(rec_num_states)
	def __call__(self, state_obj):
		return T(state_obj, self.t_amp_obj, self.rec_num_states)



class T_space_traits_cls(object):
	def __init__(self):
		self.field = np.float64
	@staticmethod
	def dot(v,w):
		return state.dot(v.t_amp_obj, w.t_amp_obj)
	@staticmethod
	def add_to(v,w,c=1):
		v.t_amp_obj.increment(w.t_amp_obj, c)
	@staticmethod
	def scale(n,v):
		v.t_amp_obj.scale(n)
	@staticmethod
	def copy(v):
		return copy.deepcopy(v)
	@staticmethod
	def act_on_vec(op,v):
		w = copy.deepcopy(v)
		w.t_amp_obj = op(v.t_amp_obj)	# This is a single tasker, only here for acting orb energy denominators
		return w
	#
	def check_member(self,v):
		pass
	def check_lin_op(self,op):
		return False
	@staticmethod
	def function_on_diags(func,op):
		raise NotImplementedError
	@staticmethod
	def back_act_on_vec(v,op):
		raise NotImplementedError
	@staticmethod
	def diagonal(op):
		raise NotImplementedError

T_space_traits = T_space_traits_cls()
