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
import numpy as np
import config
from config import configuration, transition, op_string		# so they may be imported directly from the state module without importing config



class state(object):
	def __init__(self, excite_strings, basis_vec=None):
		coefficients = np.zeros(( len(excite_strings) ))	# Numpy 1D Array
		if basis_vec is not None:  coefficients[basis_vec] = 1.
		self._data = ( deepcopy(excite_strings), coefficients )
	def __str__(self):
		excite_strings, coefficients = self._data
		string = ""
		for i in range(len(excite_strings)):  string += str(configuration(excite_strings[i],coefficients[i])) + "\n"
		return string
	def value(self):
		excite_strings, coefficients = self._data
		return (deepcopy(excite_strings), deepcopy(coefficients))
	def add_config(self, config):
		""" Add a single configuration to this state.  If config not represented in excite_strings, it will be ignored (projected off). """
		if config.is_not_null():
			excite_string,  coefficient  = config.value()
			excite_strings, coefficients = self._data
			# .index is better than searching the "list" manually because excite_strings may not be a list (though it should be indexable and iterable and have a length)
			try:                i = excite_strings.index(excite_string)
			except ValueError:  pass		# means that excite_string was not found in excite_strings
			else:               coefficients[i] += coefficient
	def increment(self, other, scaling_factor=1):
		""" Scale and add another state to this state """
		_, these_coefficients =  self._data
		_, those_coefficients = other._data
		these_coefficients += scaling_factor * those_coefficients
	def scale(self, scaling_factor):
		""" scale this state in place """
		_, coefficients =  self._data
		coefficients *= scaling_factor
	def dot(self,other):
		_, these_coefficients =  self._data
		_, those_coefficients = other._data
		return np.dot(these_coefficients,those_coefficients)
	def _op_string_action(self,ops):
		excite_strings, coefficients = self._data
		new_state = state(excite_strings)
		for pair in zip(excite_strings,coefficients):  new_state.add_config( ops | configuration(*pair) )
		return new_state

def dot(state1, state2):  return state1.dot(state2)
