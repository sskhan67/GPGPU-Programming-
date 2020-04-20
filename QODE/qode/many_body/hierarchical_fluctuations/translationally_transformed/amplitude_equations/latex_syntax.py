import inspect
from equations_transformer import *

def get_one_excitation_string(orbitals, spacing):
	upper_indices  =  '{}({})'.format( orbitals[0],   spacing[0]    )
	lower_indices  =  '{}({})'.format( orbitals[1],   spacing[1]    )
	string = 't^{{{}}}_{{{}}}'.format( upper_indices, lower_indices )
	return(string)


def get_two_excitation_string(orbitals, spacing):
	up_rgt  = '{}({})'.format( orbitals[0], spacing[0] )
	up_lft  = '{}({})'.format( orbitals[1], spacing[1] )
	low_rgt = '{}({})'.format( orbitals[2], spacing[2] )
	low_lft = '{}({})'.format( orbitals[3], spacing[3] )
	string  = 't^{{{}{}}}_{{{}{}}}'.format( up_rgt, up_lft, low_rgt, low_lft )
	return(string)


def get_restricted_sum_string(limits, index, i):
	dependent_indices = limits[i]
	dependencies = get_dependencies_string( dependent_indices )
	string       = '\sum_{{ {}{} }}'.format(index, dependencies)
	return(string)


def get_spacing_limits_string(spacing, limits):
	n_indices  = len(spacing)
	string     = ""
	if n_indices == 0:
		return string
	else:
		n_limits   = len(limits[0])
		for i in range(n_indices):
			index       = spacing[i]
			if n_limits == 2:
				if limits[i][1] == '\infty':           # This probably should be changed to 'unrestricted' string 
					lower_limit = limits[i][0]
					upper_limit = limits[i][1]
					string     += '\sum_{{{}={}}}^{{{}}}'.format( index, lower_limit, upper_limit )
				else:
					string += get_restricted_sum_string(limits, index, i)
			else:    
				string += get_restricted_sum_string(limits, index, i)
	return(string)


def get_prefactor_string(number):
	if isinstance(number, int):
		if number ==  1:
			string  = '+'
		elif number == -1:
			string  = '-'
		else:
			string  = str(number)
	else:
		integers = number.as_integer_ratio()
		string   = ""
		if integers[0] > 0:
			string  +=  "\\frac{{{}}}{{{}}}".format( integers[0], integers[1])
		else:
			string  += "-\\frac{{{}}}{{{}}}".format(-integers[0], integers[1])
	return(string)	


def get_orbitals_string(orbitals):
	if len(orbitals) > 0:
		string   = '\sum_{{{}}}'.format( ''.join(orbitals) )
	else:
		string   = ''
	return(string)

def get_hamiltonian_string(orbitals, spacing):
	n_body = len(orbitals)/2
	if n_body == 1:
		upper   = '{}({})'.format( orbitals[0], spacing[0] )
		lower   = '{}({})'.format( orbitals[1], spacing[1] )
		string  = 'f^{{{}}}_{{{}}}'.format(upper, lower)
	if n_body == 2:
		up_rgt  = '{}({})'.format( orbitals[0], spacing[0] )
		up_lft  = '{}({})'.format( orbitals[1], spacing[1] )
		low_rgt = '{}({})'.format( orbitals[2], spacing[2] )
		low_lft = '{}({})'.format( orbitals[3], spacing[3] )
		string  = 'V^{{{}{}}}_{{{}{}}}'.format(up_rgt, up_lft, low_rgt, low_lft )
	return string


def get_excitation_string(operators):
	n_operators = int( len(operators)/2 )
	string = ''
	for i in range(n_operators):
		excitation_level = int( len(operators[2*i])/2 )
		if excitation_level == 1:
			orbitals = operators[2*i]
			spacing  = operators[2*i+1]
			string  += get_one_excitation_string( orbitals, spacing )
		if excitation_level == 2:
			orbitals = operators[2*i]
			spacing  = operators[2*i+1] 
			string  += get_two_excitation_string( orbitals, spacing )
	return(string)


'''
This function takes the output from the 'normal_CC_to_TTCC' class and gives back the corresponding latex syntax of the translationally transformed coupled cluster term (TTCC).
For example:

normal_CC_term = [-1/4, [k, l, c, d], [c, i], [a, k], [d, j], [b, l], 'unrestricted' ]
TT_term        = normal_CC_to_TTCC( normal_CC_term )
print(  latex( TT_term ) )

output:

-\frac{1}{4}\sum_{o=-\infty}^{\infty}\sum_{p=-\infty}^{\infty}\sum_{q=-\infty}^{\infty}\sum_{r=-\infty}^{\infty}\sum_{ckdl}V^{k(p)l(r)}_{c(o)d(q)}t^{c(o)}_{i(n)}t^{a(m)}_{k(p)}t^{d(q)}_{j(x)}t^{b(y)}_{l(r)}

'''
def get_latex_term(TTCC_term):
	hamiltonian_orbitals = TTCC_term.hamiltonian_orbitals
	hamiltonian_spacing  = TTCC_term.hamiltonian_spacing
	prefactor            = TTCC_term.prefactor
	spacing              = TTCC_term.spacing
	limits               = TTCC_term.limits
	orbitals_sum         = TTCC_term.orbitals_sum
	ex_operators         = TTCC_term.ex_operators
	syntax               = get_prefactor_string(prefactor)
	syntax              += get_spacing_limits_string(spacing, limits)
	syntax              += get_orbitals_string(orbitals_sum)
	syntax              += get_hamiltonian_string(hamiltonian_orbitals, hamiltonian_spacing)
	syntax              += get_excitation_string(ex_operators)
	return(syntax)

def latex(terms_list):
	string = ''
	if inspect.isclass(terms_list):
		terms_list = [terms_list]
	for term in terms_list:
		term = normal_CC_to_TTCC(term)
		string += get_latex_term(term)
	print(string)	
