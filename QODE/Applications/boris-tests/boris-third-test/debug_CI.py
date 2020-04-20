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
import numpy as np
import code_CI
from class_CI import mat_dot_vec
from  code_CI import CI_mat 
from parameters_CI import states_per_oscillator, mass, k_OHconstant
from parameters_CI import k_coupling, hbar, arbitrary_mass
from numpy import linalg

H_param = {
            "mass":                  mass,
	    "states_per_oscillator": states_per_oscillator,
	    "k_OHconstant":          k_OHconstant,
	    "k_coupling":            k_coupling,
            "hbar":                  hbar,
            "arbitrary_mass":        arbitrary_mass}

# This block generates the correspoding classes to compare

# H  =  mat_dot_vec(  H_param )
# CI_matrix = CI_mat( H_param )

# Test to check matrix generated from class mat_dot_vec
'''
n_states   = states_per_oscillator ** len(mass)
I          = np.identity( n_states )
I          = np.matrix(I)
CLS_matrix = np.zeros(( n_states,n_states ))
for j in range( n_states ):
	row = I[j]
	for i in range( n_states ):
		col = I[:,i]
		CLS_matrix[i,j] = np.dot(row,H(col))

'''
#diff = np.sum( np.square( np.subtract( CI_matrix, CLS_matrix) ) )
#print("Difference between matrices is {}".format( diff ))

#print( CI_matrix  )
#print( CLS_matrix )

# Test if matrix is symmetric
#t = np.transpose( CLS_matrix )
#print(np.sum(np.square(np.subtract( CLS_matrix, t))))

# CI Matrix diagonalization and sorting eigenvalues

#sorted_eigenvals = sorted(linalg.eigvals(CI_matrix))
#print("Lowest eigenvalue from NUMPY diagonalization: {}".format(sorted_eigenvals[0]))
