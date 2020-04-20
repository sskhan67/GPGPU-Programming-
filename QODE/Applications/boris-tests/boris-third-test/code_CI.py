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
###########################################################
#                  FULL CI calculation                    #
###########################################################

#  Importing packages and functions

import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(threshold=np.nan,linewidth=1000)
import copy
from func_CI import gen_config, num_of_config, compare, CI_prod_intern, evaluate 
from numpy import linalg
from class_CI import CI_inner_product


# Filling up and printing FULL CI matrix

def CI_mat(H_param):
	n_oscillators = len(H_param["mass"])
	n_states = H_param["states_per_oscillator"] ** n_oscillators
	CI_matrix = np.zeros((n_states,n_states))
	CI_mat_element = CI_inner_product(H_param)
	for i in range(n_states):
		ket_string = gen_config(i, H_param)
		for j in range(n_states):
			bra_string = gen_config(j, H_param)
			CI_matrix[i,j] = CI_mat_element(bra_string,ket_string)
	return(CI_matrix)


