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
import inspect
from tensors import *
from letters import *
from equations_transformer import *



def get_all_free_indices(term):
	common_spacing = set( term.spacing      )  
	common_orbital = set( term.orbitals_sum )   
	full_op        = term.full_operator
	all_operators  = [full_op[j+4] for j in range(len(full_op) - 4)]
	all_orbitals   = [all_operators[2*i  ] for i in range(int(len(all_operators)/2))]
	all_free_indices = []
	for operator in all_operators:
		free_orbitals = set(operator) - common_orbital   # this still contains the common spacing indices
		free_index    = free_orbitals - common_spacing
		free_index    = list( free_index )
		if len(free_index) > 0:
			for index in free_index:   
				all_free_indices.append(index)
	return(all_free_indices)

def get_all_CC_free_indices(term):
	#common_spacing = set( term.spacing      )  
	orbitals_indices, spacing_indices = get_summation_indices(term)  # since normal CC doesn't have spacing, it will always return empty list
	#print(term, "term")
	#print(orbitals_indices, "<-------") 
	common_orbital = set( orbitals_indices )   
	full_op        = term 
	all_operators  = [full_op[j+1] for j in range(len(full_op) - 2)]
	#all_orbitals   = [all_operators[2*i  ] for i in range(int(len(all_operators)/2))]
	all_free_indices = []
	for operator in all_operators:
		free_orbitals = set(operator) - common_orbital   # this still contains the common spacing indices
		#free_index    = free_orbitals - common_spacing
		free_index    = list( free_orbitals )
		if len(free_index) > 0:
			for index in free_index:   
				all_free_indices.append(index)
	return(all_free_indices)

def get_free_spacings(term):
	common_spacing = set( term.spacing )
	full_op        = term.full_operator
	all_operators  = [full_op[j+4] for j in range(len(full_op) - 4)]
	all_spacings   = [all_operators[2*i+1] for i in range(int(len(all_operators)/2))]
	free_spacings  = []
	for operator in all_spacings:
		free_oper_spacings     = set(operator) - common_spacing
		free_oper_spacings     = list(free_oper_spacings)
		if len(free_oper_spacings) > 0:
			free_spacings += free_oper_spacings
	return(free_spacings)




def get_orbital_position(excit_operators, orbital):     
	for operator in excit_operators:
		if orbital in operator:
			index    = operator.index(orbital)
			position = index*2
			if len(operator) == 2:
				string   = 't_ai'
			if len(operator) == 4:
				string   = 't_abij'
			for i in range(position):
				string += '[0]'
			return(string)


def get_orbital_identity(orbital):
	virtuals = ['a', 'b', 'c', 'd', 'e']
	occupied = ['i', 'j', 'k', 'l']
	if orbital in virtuals:
		output_str = 'n_virtuals'
	if orbital in occupied:
		output_str = 'n_occupied'
	return output_str


def get_imports_and_intervals(interval, transformed_terms):
	if isinstance(transformed_terms[0], TTCC_term):  
		output_str  = '\nfrom tensors import *\n\n\n'
		output_str += 'L = {}\n'.format(interval)
		output_str += get_function_heading(transformed_terms[0])
		output_str += '\n    m=0\n'
		#output_str += 'result = 0\n\n\n'
	else:
		output_str  = get_CC_function_heading(transformed_terms[0])
	return output_str 


def get_file_heading(transformed_terms):
	if isinstance(transformed_terms[0], TTCC_term):  
		output_str  = '\n\n'
		#output_str += 'L = {}\n'.format(interval)
		output_str += get_function_heading(transformed_terms[0])
		output_str += '\n    m=0\n'
		#output_str += 'result = 0\n\n\n'
	else:
		output_str  = '\n\n'
		output_str += get_CC_function_heading(transformed_terms[0])
	return output_str 

# this function does not use ex_operators anymore, you need to change it

# this function does not use ex_operators anymore, you need to change it
def get_orbital_space_identity(orbitals_sum):  # you may want to number the limits by 1, instead of 2
	limits         = []
	#limit_orbitals = 1
	for orbital in orbitals_sum:    
		#orbital_position = get_orbital_position(ex_operators, orbital) # position in terms of '[0][0]'
		orbital_identity = get_orbital_identity(orbital)  # in terms of virtual/occupied
		#tuple            = (limit_orbitals, orbital_identity)
		limits.append(orbital_identity) 
		#limit_orbitals  += 2
	return(limits)


def get_upper_limits_string(orbitals_sum, ex_operators): # probably this functions is now used, consider deletion
	limits     = get_orbital_domains(orbitals_sum, ex_operators)
	n_limits   = len(limits)
	output_str = ''
	for i in range(n_limits):
		output_str += '    upper_lim{} =  {} \n'.format( limits[i][0], limits[i][1] )
	return(output_str)


def get_initial_values_string(term):
	output_str          = ''
	all_free_indices    = get_all_free_indices(term)
	print(all_free_indices, 'free indices')
	for index in all_free_indices:
		output_str += '    {} = 0\n'.format(index)
	return(output_str)

def get_function_heading(term):
	#print(term.full_operator)
	#output_str        = 'def delta_t('
	all_free_indices  = get_all_free_indices(term)
	excitation_level = int(len(all_free_indices)/4)
	if excitation_level == 1 :
		output_str = 'def delta_t(a,n,i,m,'
	if excitation_level == 2 :
		output_str = 'def delta_t(i,j,a,b,n_1,n_2,n_3,'
	if excitation_level == 0 :
		output_str = "def delta_E("
	#for index in all_free_indices:
	#	output_str += '{},'.format(index)
	output_str += 'L, n_virtuals, n_occupied, f_pq, v_pqrs, t_ai, t_abij):\n'
	output_str += '    result = 0\n'
	return output_str


def get_CC_function_heading(term):
	all_free_indices  = get_all_CC_free_indices(term)
	excitation_level = int(len(all_free_indices)/2)
	if excitation_level == 1 :
		output_str = 'def delta_t(a, i,'
	if excitation_level == 2 :
		output_str = 'def delta_t(i, j, a, b,'
	if excitation_level == 0 :
		output_str = "def delta_E("
	#for index in all_free_indices:
	#	output_str += '{},'.format(index)
	output_str += 'n_virtuals, n_occupied, f_pq, v_pqrs, t_ai, t_abij):\n'
	output_str += '    result = 0\n'
	return output_str


'''
def get_intervals(term, white_space, past_index):
	string = ''
	dependencies = term[2]
	n_dependecies = len(term[2])
        # all_free_indices, free_spacings = get_free_indices(term) 
	for operator in dependencies:
		for index in operator:
			if index not in free_spacings:
				string += '{}L_{} = [i for i in range({}-L, {}+L+1)]\n'.format(white_space, center, center, center)
	return(string)
'''

def get_intersections(current_spacing, dependencies, white_space):
	l = 0
	for center in dependencies:
		if l == 0:
			string  = '    {}domain_{} = list( set(L_{}) '.format(white_space, current_spacing, center)
		else:
			string += ' & set(L_{}) '.format(center)
		l += 1
	string += ')\n'
	return(string)

def get_orbital_loop(white_space2, orbitals_sum, orbital_domains, i):
	if orbital_domains[i] == 'n_virtuals':
		output_str = '    {}for {} in range(n_occupied,{}):\n'.format(white_space2,  orbitals_sum[i], orbital_domains[i])
	else:
		output_str = '    {}for {} in range({}):\n'.format(white_space2,  orbitals_sum[i], orbital_domains[i])
	return output_str 
	
def get_first_loop(white_space, white_space2, spacing, orbitals_sum, orbital_domains, i):
	#orbital_identity = get_orbital_identity(orbital)  # in terms of virtual/occupied
	output_str  = ''
	output_str += '    {}for {} in domain_{}:\n'.format(         white_space,   spacing[i],      spacing[i])
	output_str += '    {}for {} in range({}):\n'.format(white_space2,  orbitals_sum[i], orbital_domains[i])
	#output_str +=  get_orbital_loop(white_space2, orbitals_sum, orbital_domains, i)
	return(output_str)

def get_first_CC_loop(white_space, white_space2, orbitals_sum, orbital_domains, i):
	#orbital_identity = get_orbital_identity(orbital)  # in terms of virtual/occupied
	output_str  = ''
	#output_str += '    {}for {} in domain_{}:\n'.format(         white_space,   spacing[i],      spacing[i])
	output_str += '{}for {} in range({}):\n'.format(white_space2,  orbitals_sum[i], orbital_domains[i])
	#output_str +=  get_orbital_loop(white_space2, orbitals_sum, orbital_domains, i)
	return(output_str)


def get_subsequent_loops(white_space, white_space2, spacing, orbitals_sum, orbital_domains, past_index, limits, i):
	output_str  = ''
	output_str += '    {}L_{} = [i for i in range({}-L, {}+L+1)]\n'.format(white_space, past_index, past_index, past_index)
	output_str += get_intersections( spacing[i], limits[i], white_space)
	output_str += '    {}for {} in domain_{}:\n'.format(         white_space,   spacing[i],      spacing[i])
	output_str += '    {}for {} in range({}):\n'.format(white_space2,  orbitals_sum[i], orbital_domains[i])
	#output_str += get_orbital_loop(white_space2, orbitals_sum, orbital_domains, i)
	return(output_str)

def get_subsequent_CC_loops(white_space, white_space2, orbitals_sum, orbital_domains, i):
	output_str  = ''
	#output_str += '    {}L_{} = [i for i in range({}-L, {}+L+1)]\n'.format(white_space, past_index, past_index, past_index)
	#output_str += get_intersections( spacing[i], limits[i], white_space)
	#output_str += '    {}for {} in domain_{}:\n'.format(         white_space,   spacing[i],      spacing[i])
	output_str += '{}for {} in range({}):\n'.format(white_space2,  orbitals_sum[i], orbital_domains[i])
	#output_str += get_orbital_loop(white_space2, orbitals_sum, orbital_domains, i)
	return(output_str)


def get_loops_string(spacing, limits, orbitals_sum):
	output_str   = ''
	n_loops      = len( orbitals_sum )
	counter      = 0
	white_space3 = ''
	orbital_domains = get_orbital_space_identity( orbitals_sum )
	for i in range(n_loops):
		white_space = counter    *'    '
		white_space2= (counter+1)*'    '
		white_space3= (counter+2)*'    '
		lim_number  = counter + 1
		if counter == 0:
			output_str += get_first_loop(white_space, white_space2, spacing, orbitals_sum, orbital_domains, i)
			counter    += 2
		else:
			past_index  = spacing[i-1]
			output_str += get_subsequent_loops(white_space, white_space2, spacing, orbitals_sum, orbital_domains, past_index, limits, i)
			counter    += 2
	return(output_str, white_space3)

def get_CC_loops_string(orbitals_sum):
	output_str   = ''
	n_loops      = len( orbitals_sum )
	counter      = 0
	white_space3 = ''
	orbital_domains = get_orbital_space_identity( orbitals_sum )
	for i in range(n_loops):
		white_space = counter    *'    '
		white_space2= (counter+1)*'    '
		white_space3= (counter+2)*'    '
		lim_number  = counter + 1
		if counter == 0:
			output_str += get_first_CC_loop(white_space, white_space2, orbitals_sum, orbital_domains, i)
			counter    += 2
		else:
			#past_index  = spacing[i-1]
			output_str += get_subsequent_CC_loops(white_space, white_space2, orbitals_sum, orbital_domains, i)
			counter    += 2
	return(output_str, white_space3)


def get_orbital_realignment(orbitals_sum, white_space):
	virtual_orbitals = ['a', 'b', 'c', 'd']
	output_str = ''
	for orbital in virtual_orbitals:
		if orbital in orbitals_sum:
			output_str += '{}    {} = {} + n_occupied\n'.format(white_space, orbital, orbital)
	return output_str

def add_comma_or_parenthesis(index, limit):
	if index < limit:
		output_str = ', '
	else:
		output_str = ')'
	return(output_str)

'''
				if l == 0:
					output_str += ', '
					l += 1
				else:
					output_str += ')'
'''

def get_hamiltonian_product_string(prefactor, H_orbitals, H_spacings, white_space):
	virtuals = ['a', 'b', 'c', 'd']
	output_str = '    {}result += {}*'.format(white_space, prefactor)
	if len(H_orbitals) == 2:
		output_str += "f_pq("
		l = 0
		parenthesis_position = 1
		for i in range(2):
			if H_orbitals[i] in virtuals:
				output_str += '{} + n_occupied, {}'.format(H_orbitals[i], H_spacings[i])
				output_str += add_comma_or_parenthesis(l, parenthesis_position)  
			else:
				output_str += '{}, {}'.format(H_orbitals[i], H_spacings[i])
				output_str += add_comma_or_parenthesis(l, parenthesis_position)
			l += 1
	else:
		output_str += "v_pqrs("
		m = 0
		parenthesis_position = 3
		for i in range(4):
			if H_orbitals[i] in virtuals:
				output_str += '{} + n_occupied, {}'.format(H_orbitals[i], H_spacings[i])
				output_str += add_comma_or_parenthesis(m, parenthesis_position)
			else:
				output_str += '{}, {}'.format(H_orbitals[i], H_spacings[i])
				output_str += add_comma_or_parenthesis(m, parenthesis_position)
			m += 1
	return(output_str)

def get_CC_hamiltonian_product_string(prefactor, H_orbitals, white_space):
	virtuals = ['a', 'b', 'c', 'd']
	output_str = '    {}result += {}*'.format(white_space, prefactor)
	if len(H_orbitals) == 2:
		output_str += "f_pq("
		l = 0
		parenthesis_position = 1
		for i in range(2):
			if H_orbitals[i] in virtuals:
				output_str += '{} + n_occupied'.format(H_orbitals[i])
				output_str += add_comma_or_parenthesis(l, parenthesis_position)  
			else:
				output_str += '{}'.format(H_orbitals[i])
				output_str += add_comma_or_parenthesis(l, parenthesis_position)
			l += 1
	else:
		output_str += "v_pqrs("
		m = 0
		parenthesis_position = 3
		for i in range(4):
			if H_orbitals[i] in virtuals:
				output_str += '{} + n_occupied'.format(H_orbitals[i])
				output_str += add_comma_or_parenthesis(m, parenthesis_position)
			else:
				output_str += '{}'.format(H_orbitals[i])
				output_str += add_comma_or_parenthesis(m, parenthesis_position)
			m += 1
	return(output_str)


'''
	output_str = ''
	if len(H_orbitals) == 2:
		T = ( white_space, H_orbitals[0], H_spacings[0], H_orbitals[1] , H_spacings[1], prefactor  )
		output_str += '    {}result += {}*f_pq({}, {}, {}, {})'.format(T[0], T[-1], T[1], T[2], T[3], T[4])
	else:
		T = ( white_space, H_orbitals[0], H_spacings[0], H_orbitals[1], H_spacings[1], H_orbitals[2], H_spacings[2], H_orbitals[3], H_spacings[3], prefactor)
		output_str += '    {}result += {}*v_pqrs({}, {}, {}, {}, {}, {}, {}, {})'.format(T[0], T[-1], T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8])
	return(output_str)
'''

def get_exc_indices_pair(n_index, excit_operators, i, j, counter):
	output_str = ''
	if counter == 0:
		if n_index == 2:
			output_str += '*t_ai({}, {}'.format(excit_operators[2*i][j], excit_operators[2*i+1][j])
		if n_index == 4:
			output_str += '*t_abij({}, {} '.format(excit_operators[2*i][j], excit_operators[2*i+1][j])
	else:
		output_str += ', {}, {}'.format(excit_operators[2*i][j], excit_operators[2*i+1][j])
	return(output_str)


def get_CC_exc_indices_pair(n_index, excit_operators, i, j, counter):
	output_str = ''
	if counter == 0:
		if n_index == 2:
			output_str += '*t_ai({}'.format(excit_operators[i][j])
		if n_index == 4:
			output_str += '*t_abij({}'.format(excit_operators[i][j])
	else:
		output_str += ', {}'.format(excit_operators[i][j])
	return(output_str)



def get_excitation_product_string(excit_operators):
	output_str = ''
	n_pairs    = int(len(excit_operators)/2)     # number of pairs of operators
	for i in range(n_pairs):
		n_index = len(excit_operators[2*i])
		counter = 0
		for j in range(n_index):
			output_str += get_exc_indices_pair(n_index, excit_operators, i, j, counter)
			counter    += 1
		output_str += ')'
	return(output_str)

def get_CC_excitation_product_string(excit_operators):
	output_str = ''
	n_operators   = len(excit_operators)     # number of operators
	for i in range(n_operators):
		n_index = len(excit_operators[i])
		counter = 0
		for j in range(n_index):
			output_str += get_CC_exc_indices_pair(n_index, excit_operators, i, j, counter)
			counter    += 1
		output_str += ')'
	return(output_str)

def get_free_intervals_string(term):
	free_spacings = get_free_spacings(term)
	#print('free spacings', free_spacings)
	string        = ''
	n_free        = len(free_spacings)
	if n_free > 0:
		for index in free_spacings:
			string += '    L_{} = [i for i in range({}-L, {}+L+1)]\n'.format(index, index, index)
	#print('term.limits=', term.limits )
	if  term.limits:                   # check if list is filled up
		if not term.limits[0]:     # check if list is empy
			string += '    L_0 = [i for i in range(-L, L+1)]\n'
	return(string)

def get_initial_intersections_string(spacings, limits):
	string         = ''
	n_dependencies = len(limits)
	if n_dependencies == 0:
		return string
	initial_dependencies = limits[0]
	n_indices            = len(initial_dependencies)
	initial_spacing      = spacings[0]
	if n_indices  == 0:
		center = 0
		string = '    domain_{} = L_0'.format(initial_spacing)
	else: 
		if initial_dependencies[0] == '-\infty':
			print(  "IDIOT, you cannot program an infinite term" )
		if n_indices == 1:
			center = initial_dependencies[0]
			string = '    domain_{} = L_{}'.format(initial_spacing, center)
		else:
			no_space = ''                     					       
			string  += get_intersections(initial_spacing, initial_dependencies, no_space)  
	return(string)


def get_TTCC_summations(TTCC_term):     
	prefactor                = TTCC_term.prefactor
	spacing                  = TTCC_term.spacing
	limits                   = TTCC_term.limits
	orbitals_sum             = TTCC_term.orbitals_sum
	hamiltonian_orbitals     = TTCC_term.hamiltonian_orbitals
	hamiltonian_spacing      = TTCC_term.hamiltonian_spacing
	ex_operators             = TTCC_term.ex_operators
	loop_string, white_space = get_loops_string(spacing, limits, orbitals_sum)
	#output_str               = get_initial_values_string(TTCC_term)
	output_str               = '\n\n'
	#output_str              += get_upper_limits_string(orbitals_sum, ex_operators)
	output_str              += get_free_intervals_string(TTCC_term)
	output_str              += get_initial_intersections_string(spacing, limits)
	output_str              += '\n\n'
	output_str              += loop_string
	#output_str              += get_orbital_realignment(orbitals_sum, white_space)
	output_str              += get_hamiltonian_product_string(prefactor, hamiltonian_orbitals, hamiltonian_spacing, white_space)
	output_str              += get_excitation_product_string(ex_operators)
	output_str              += '\n\n\n'	
	return(output_str)

def get_CC_summations(CC_term):     
	prefactor                = CC_term[0]
	orbitals_sum,spacing_sum = get_summation_indices(CC_term)  # since normal CC doesn't have spacing, it will return an empty list
	hamiltonian_orbitals     = CC_term[1]
	ex_operators             = CC_term[2:len(CC_term)-1]  
	loop_string, white_space = get_CC_loops_string(orbitals_sum)
	output_str               = '\n\n'
	#output_str              += get_upper_limits_string(orbitals_sum, ex_operators)
	#output_str              += get_free_intervals_string(TTCC_term)
	#output_str              += get_initial_intersections_string(spacing, limits)
	output_str              += '\n\n'
	output_str              += loop_string
	#output_str              += get_orbital_realignment(orbitals_sum, white_space)
	output_str              += get_CC_hamiltonian_product_string(prefactor, hamiltonian_orbitals, white_space)
	output_str              += get_CC_excitation_product_string(ex_operators)
	output_str              += '\n\n\n'	
	return(output_str)

def get_summations_code(term):
	if isinstance(term, TTCC_term):
		code_string = get_TTCC_summations(term)
	else:
		code_string = get_CC_summations(term)
	return code_string


def code(term_list, file_name, CC_mode):  
	#file_name = "generated_CC_equation_code.py" 
	if inspect.isclass(term_list):
		term_list = [term_list]
	f = open(file_name, 'a')
	transformed_terms = []
	#print(term_list, "term list")
	for term in term_list:
		if CC_mode == 'normal CC':
			pass
		elif CC_mode == 'translationally transformed CC': # switch from normal CC to TTCC
			term = normal_CC_to_TTCC(term)  
		else:
			print('you did not specify CC mode, so it will be normal CC as default')
		transformed_terms.append(term)
	#f.write(get_imports_and_intervals(interval_limit, transformed_terms))
	f.write(get_file_heading(transformed_terms))
	for trans_term in transformed_terms:
		f.write(get_summations_code(trans_term))
	f.write('\n    return(result)\n')
	f.close()
	print('"{}" has been generated'.format(file_name))



'''
# Test of isolated module
f = open('generated_code.py', 'a')
#normal_CC_term = [-1/4, [k, l, c, d], [c, i], [a, k], [d, j], [b, l], 'restricted' ]
normal_CC_term = [+1/2, [i, j, a, b], [a, i,], [b, j], 'restricted']
term        = normal_CC_to_TTCC( normal_CC_term )
f.write( code(term) )
f.close()
'''
