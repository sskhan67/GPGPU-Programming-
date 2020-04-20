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
from copy import copy
from qode.fermion_field.occ_strings import all_occ_strings as excited_molecules		# a coincident usage of same logic


# The recursive workhorse and its user-friendly wrapper

def _make_excite_strings(excite_string,n_states):
	"""\
	This is a workhorse function, not likely to be of use to a caller outside of this file.
	Given a configuration (list of states) for N molecules, it returns a list of configurations for
	N+1 molecules, where the additional molecule can be in any one of n_states[0] states.
	"""
	the_list = []
	if len(n_states)==0:
		return [excite_string]
	else:
		for n in range(n_states[0]):
			excite_string_temp = copy(excite_string)
			excite_string_temp += [n]
			the_list += _make_excite_strings(excite_string_temp,n_states[1:])
	return the_list

def all_excite_strings(numbers_of_states):
	"""\
	This wrapper is suitable for generating all FCI occupation strings for len(numbers_of_states) molecules,
	where each molecule has an associated number of states given by the corresponding element of numbers_of_states.
	A list of lists is returned, where the internal lists give the state of each molecule in a given configuration.
	"""
	return _make_excite_strings([],numbers_of_states)



# Helpers for generating occupation strings of a certain excitation order, given a reference implied by the occupied set

def _n_tuple_substitutions(n, numbers_of_states):
	"""\
	Generates all occupation strings associated with given exciation order n above a reference,		# written quickly and never tested!!
	taken to be all molecules in state 0.  Returns the result as a list of occupation strings,
	each of which is a list of the state for each molecule.
	"""
	all_configs = []
	excitation_groups = excited_molecules(len(numbers_of_states), n)
	for excitation_group in excitation_groups:
		sub_list = []
		for M,molecule in enumerate(excitation_group):
			if molecule:  sublist += numbers_of_states[M]
		excite_strings = all_excite_strings(sub_list)
		for excite_string in excite_strings:
			i = 0
			config = []
			for M,molecule in enumerate(excitation_group):
				if molecule:
					config += [excite_string[i]]
					i += 1
				else:
					config += [0]
			all_configs += [config]
	return all_configs



# CI occupation string sets

def CI(orders, numbers_of_states):
	"""\
	occ_orbs and vrt_orbs are lists containing the indices of occupied and virtual spin orbitals in any order.

	This returns a list of occupation strings (each of which is a list of boolean occupancies) representing all excitations
	of all orders in orders (a list of integers) above a reference implied by the occupied set, in an orbital space that is
	the union of the occupied and virtual sets.

	The blocking by excitation order will be in the order given in the input list, but the ordering of occupation
	strings associated with a given excitation order is defined internally by the algorithm.

	This basically functions as a more sophisticated wrapper for _n_tuple_substitutions().
	"""
	excite_strings = []
	for n in orders:
		if n>len(numbers_of_states):  raise Exception("CI excitation order exceeds available number of molecules")
		excite_strings += _n_tuple_substitutions(n, numbers_of_states)
	return excite_strings

def CIS(   numbers_of_states):  return CI([0,1],      numbers_of_states)
def CISD(  numbers_of_states):  return CI([0,1,2],    numbers_of_states)
def CISDT( numbers_of_states):  return CI([0,1,2,3],  numbers_of_states)
def CISDTQ(numbers_of_states):  return CI([0,1,2,3,4],numbers_of_states)

def FCI(numbers_of_states):	# Not blocked by excitation order, and pays no heed to orbital ordering!
	return all_excite_strings(numbers_of_states)
