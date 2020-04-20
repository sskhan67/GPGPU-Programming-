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
import qode
import numpy
from debug_CI import H_param
from my_space import space_traits
from class_CI import mat_dot_vec
from code_CI  import CI_mat
from numpy import linalg

def energy_to_file( eigenval, oscillators, states, coupling):
    output_str  = ""
    output_str += "Oscillators: {}".format(oscillators)
    output_str += "  States:  {}".format(states)
    output_str += "  Coupling: {}".format(coupling)
    output_str += "  Lowest eigenval Numpy diagonal:    {}\n".format(eigenval) 
    return(output_str) 

for j in range(3,7):
	for i in range(3,7):
	  for k in [0]:
		  H_param["states_per_oscillator"] = i
		  number_of_oscillators            = j
		  H_param["mass"]         = number_of_oscillators * [ H_param["mass"][0] ]
		  H_param["k_OHconstant"] = number_of_oscillators * [ H_param["k_OHconstant"][0] ]
		  H_param["k_coupling"]   = number_of_oscillators * [ k ]
		  CI_matrix = CI_mat( H_param )
		  sorted_eigenvals = sorted(linalg.eigvals(CI_matrix))
		  f = open('results_numpy_2.out', 'a')
		  f.write( energy_to_file(sorted_eigenvals[0], j, i, k) )
		  f.close()
