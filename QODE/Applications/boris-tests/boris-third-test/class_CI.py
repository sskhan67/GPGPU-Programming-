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
########################################################
#                       Class module                   #
########################################################
import sys
import numpy
#numpy.set_printoptions(threshold=sys.maxsize)
numpy.set_printoptions(threshold=numpy.nan,linewidth=1000)
import copy
from numpy import linalg
from func_CI import eval_up_up, eval_dwn_dwn, eval_dwn_up, eval_up_dwn, eval_diagonal
from func_CI import gen_config, compare, evaluate, num_of_config, CI_prod_intern


class CI_inner_product(object):
	def __init__(self, H_param):
		self.H_param        = H_param
	def __call__(self, bra, ket):
		mat_element = CI_prod_intern( bra, ket, self.H_param)
		return(mat_element) 



# This class generates the different non-zero elements using ONLY ONE "for loop"
# and uses MY FUNCTION "evaluate".
#
# Details:
# This class uses MY FUNCTION "evaluate" for redundant lines of code in evaluating 
# each new generated state.
# The FUNCTION "evaluate" does all the arimetic operations each time a new state 
# is generated.
#
class mat_dot_vec(object):
	def __init__(self, H_param):
		self.H_param       = H_param
		self.n_oscillators = len(self.H_param["mass"])
		self.n_states      = H_param["states_per_oscillator"] ** self.n_oscillators
	def __call__(self,vec):
		CI_mat_el =  CI_inner_product(self.H_param)
		final_vec =  numpy.zeros((len(vec),1))
		HO_lim    =  len(self.H_param["mass"])       - 1
		up_lim    =  self.H_param["states_per_oscillator"] - 1
		for i in range(self.n_states):
			ket_str = gen_config(i, self.H_param)
			eval_diagonal( up_lim, final_vec, i, vec, ket_str, CI_mat_el, self.H_param)
			for j in range(-1, HO_lim):
				eval_up_up(  up_lim, final_vec, i, vec, ket_str, j, CI_mat_el, self.H_param)
				eval_dwn_dwn(up_lim, final_vec, i, vec, ket_str, j, CI_mat_el, self.H_param)
				eval_dwn_up( up_lim, final_vec, i, vec, ket_str, j, CI_mat_el, self.H_param)
				eval_up_dwn( up_lim, final_vec, i, vec, ket_str, j, CI_mat_el, self.H_param) 
		return(final_vec)

