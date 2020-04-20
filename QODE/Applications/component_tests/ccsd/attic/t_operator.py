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
from copy import deepcopy
import numpy
from qode.fermion_field import state





def T1(state_obj, t_amp_obj, occ_orbs, vrt_orbs, offset):
	# Take any state_obj but make a new_state
	#
	# T1 = \sum_{a i} t_ia a+ i
	# 
	#   t_ia are the coeffs in state_obj but only the SINGLE excitation amplitudes.
	#
	num_occ_orb = len(occ_orbs)
	num_vrt_orb = len(vrt_orbs)
	occ_strings,_ = state_obj.value()
	new_state = state.state(occ_strings)
	_,t_amps  = t_amp_obj.value()
	#
	#
	single_idx = offset
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		#
		#
		for a in range(num_vrt_orb):
			orb_a = vrt_orbs[a]
			#
			#
			op = state.op_string( state.create(orb_a), state.annihilate(orb_i) )
			new_state.increment( op(state_obj), t_amps[single_idx] )
			single_idx += 1
	#
	#
	return new_state



	

def T2(state_obj, t_amp_obj, occ_orbs, vrt_orbs, offset):
	# Take any state_obj but make a new_state
	#
	# T2 = \sum_{a<b, i<j} t_{ijab} a+ i b+ j 
	#
	#   t_{ijab} are the coeffs in state_obj but only the DOUBLE excitation amplitudes.
	#
	num_occ_orb = len(occ_orbs)
	num_vrt_orb = len(vrt_orbs)
	occ_strings,_ = state_obj.value()
	new_state = state.state(occ_strings)
	_,t_amps  = t_amp_obj.value()
	#
	#
	#
	double_idx = offset
	for i in range(num_occ_orb):
		orb_i = occ_orbs[i]
		for j in range(i+1, num_occ_orb):
			orb_j = occ_orbs[j]
			#
			#
			for a in range(num_vrt_orb):
				orb_a = vrt_orbs[a]
				for b in range(a+1, num_vrt_orb):
					orb_b = vrt_orbs[b]
					#
					#
					op = state.op_string( state.create(orb_a), state.annihilate(orb_i), state.create(orb_b), state.annihilate(orb_j) )
					new_state.increment( op(state_obj), t_amps[double_idx] )
					double_idx += 1
	#
	#
	return new_state




def T_with_complicated_arguments(state_obj, trivial_list_of_idx, t_amp_obj, occ_orbs, vrt_orbs, offsets, trivial_h_mat, trivial_V_mat):
	# Take any state_obj but make a new_state
	#
	# Call T1() and T2()
	#  T|CC> = T1|CC> + T2|CC>
	#
	num_R, num_RS = offsets
	new_state =         T1(state_obj, t_amp_obj, occ_orbs, vrt_orbs, num_R)
	new_state.increment(T2(state_obj, t_amp_obj, occ_orbs, vrt_orbs, num_RS))
	return new_state



# A very thin wrapper that mimics a functional programming style nested function definition, because multiprocessing chokes on nested function defs.
class T_operator(object):
	def __init__(self, t_amp_obj, occ_orbs, vrt_orbs, offsets):
		self.t_amp_obj = t_amp_obj
		self.occ_orbs  = occ_orbs
		self.vrt_orbs  = vrt_orbs
		self.offsets   = offsets
	def __call__(self,state_obj):
		return T_with_complicated_arguments(state_obj, None, self.t_amp_obj, self.occ_orbs, self.vrt_orbs, self.offsets, None, None)



class T_space_traits_class(object):
	def __init__(self):
		self.field = numpy.float64
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
		return deepcopy(v)
	@staticmethod
	def act_on_vec(op,v):
		w = deepcopy(v)
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

T_space_traits = T_space_traits_class()
