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
from letters           import *
from index_restriction import get_restrictions
from index_restriction import get_dependencies_string



class TTCC_term:    # This class will be where all information will be extrated from.
	def __init__(self, parameters):
		self.prefactor             = parameters[0] 
		self.spacing               = parameters[1]
		self.limits                = parameters[2]
		self.orbitals_sum          = parameters[3]
		self.hamiltonian_orbitals  = parameters[4]
		self.hamiltonian_spacing   = parameters[5]
		self.ex_operators          = parameters[6:]
		self.full_operator         = parameters


def get_summation_indices(CC_term):
	operators                = [CC_term[j+1] for j in range(len(CC_term) - 2) ]  # A dirty way to discard the prefactor and the 'restriction' string, but it works. 
	possible_spacing_indices = [o, p, q, r, s, t, u]
	common_orbitals          = []
	hamiltonian              = set(operators[0])
	n                        = len(operators)
	for i in range( n-1 ):
		common_orbitals += list( hamiltonian.intersection(operators[i+1]))
	N = len(common_orbitals)
	common_spacings = possible_spacing_indices[:N]
	return(common_orbitals, common_spacings)


def get_spacing_limits(CC_term, spacing_sum, all_spacings):
	n_limits        = len(spacing_sum)
	if CC_term[-1] == 'unrestricted':
		possible_limits = [[moo, oo],[moo, oo], [moo, oo], [moo, oo]]    # moo = -oo = minus infinity
		limits          = possible_limits[:n_limits]
	if CC_term[-1] == 'restricted':
		limits          = get_restrictions(spacing_sum, all_spacings)            
	return limits


def get_hamiltonian_spacings( CC_term, 
			      common_orbitals, 
			      common_spacings, 
			      possible_free_indices
			      ):
	hamiltonian_orbitals = CC_term[1]
	spacings             = get_spacing_indices( hamiltonian_orbitals, 
						    common_spacings, 
						    common_orbitals, 
						    possible_free_indices )
	return spacings


def get_spacing_index( orbital, 
		       oper_common_orbitals, 
		       free_indices, 
		       possible_free_indices, 
		       common_spacings, 
		       common_orbitals
		       ):
	if orbital in oper_common_orbitals:   # oper_common_orbitals = operator's common orbitals
		position = common_orbitals.index(orbital)
		index    = common_spacings[position]
	if orbital in free_indices:
		#index = possible_free_indices.pop()  # next modification/assign directly indices
		if   orbital == 'a':
			index = 'n_2'
		elif orbital == 'b':
			index = 'n_3'
		elif orbital == 'j':
			index = 'n_1'
		elif orbital == 'i':
			index = 'm'  
		#print(orbital, '<-- orbital') 
	return(index)
		
	
def get_spacing_indices( operator_orbitals, 
			 common_spacings, 
			 common_orbitals, 
			 possible_free_indices
			 ): 
	oper_common_orbitals = set(operator_orbitals) & set(common_orbitals)  # & = intersection between sets
	free_indices         = set(operator_orbitals) - set(common_orbitals)
	spacing_indices      = []
	for orbital in operator_orbitals:
		index = get_spacing_index( orbital, 
					   oper_common_orbitals, 
					   free_indices, 
					   possible_free_indices, 
					   common_spacings, 
					   common_orbitals )     # for a given operator orbital, it assigns a common or free index 
		spacing_indices.append(index)
	return spacing_indices


def get_exc_operators_spacings( CC_term, 
				common_spacings, 
				common_orbitals, 
				possible_free_indices ):
	exc_operators = [CC_term[j+2] for j in range(len(CC_term) - 3) ]  # A dirty way to get only the excitation operators, but it works. 
	result        = []
	spacings      = []
	for operator in exc_operators:
		exc_indices = get_spacing_indices( operator, 
						   common_spacings, 
						   common_orbitals, 
						   possible_free_indices)
		result.append(operator)
		result.append(exc_indices)
		spacings.append(exc_indices)
	return result, spacings



'''
Things to keep in mind about the class input:

 * This function only accepts its inputs as a list of lists. Although, the first and last arguments are not contained in a list.
 * The first argument is the prefactor, so it has to be a number. 
 * The last argument is one of two strings('unrestricted' or 'restricted'). It controls the restrictions on the spacing limits: 'unrestricted' == -oo to +oo; 'restricted' == restrictions given by get_restrictions().
 * No need to input the summation indices of orbitals. The code will deduce them automatically.
 * The indices of a single integral or amplitude are contained in a list.

Example of input 
for the term in conventional coupled cluster(latex syntax):

\frac{1}{4}\sum_{klcd}V^{kl}_{cd}t^{c}_{i}t^{a}_{k}t^{d}_{n}t^{b}_{l}

the input syntax is:

normal_term = [1/4, [k, l, c, d], [c, i], [a, k], [d, j], [b, l], 'unrestricted' ]

'''
def normal_CC_to_TTCC( normal_CC_term ):           # TTCC = Translationally Transformed Coupled Cluster term
	#possible_free_indices    = [n, m, w, x, y, z] # Used indices will be popped out from the list
	possible_free_indices    = [n_1, m, n_2, n_3, n_4, n_5] # Used indices will be popped out from the list
	possible_free_indices.reverse()            # There is no direct way to delete elements from the left, hence the reverse()
	#print(normal_CC_term, "normal_CC_term")
	prefactor                    = normal_CC_term[0]
	orbital_sum, spacing_sum     = get_summation_indices(normal_CC_term)
	hamiltonian_orbitals         = normal_CC_term[1]
	hamiltonian_spacings         = get_hamiltonian_spacings(normal_CC_term, orbital_sum, spacing_sum, possible_free_indices)
	exc_operators, exc_spacings  = get_exc_operators_spacings(normal_CC_term, spacing_sum, orbital_sum, possible_free_indices)
	all_spacings                 = [hamiltonian_spacings] + exc_spacings
	spacing_limits               = get_spacing_limits(normal_CC_term, spacing_sum, all_spacings)
	TTCC_term_inputs             = [prefactor, spacing_sum, spacing_limits, orbital_sum, hamiltonian_orbitals, hamiltonian_spacings]
	TTCC_term_inputs            += exc_operators
	#print(TTCC_term_inputs, 'Transformed terms')
	Translationally_Transformed_term = TTCC_term( TTCC_term_inputs )     # Class to wrap and label the TTCC inputs
	return( Translationally_Transformed_term )

