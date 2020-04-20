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
import copy



class configuration(object):
	"""\
	The main job of this class is to store an excitation string and be able to resolve the action of a
	transition operator on that string.  We allow a coefficient to be passed along, since that will also be a
	convenience for functions that use this.

	In the interest of being explicit, this class expects you to test whether or not it represents the Null
	configuration from the outside, because it will raise an Exception if you try to look inside Null.
	"""
	def __init__(self, excite_string, coefficient=1):
		# excite_string = [M,N,...],  ... a list of the states of individual molecules
		self._excite_string  = excite_string
		self._coefficient    = coefficient
	def __str__(self):
		string = "{:+6.3f}   x   | ".format(self._coefficient)		# do for Null too (coeff=0), as cheap way of keeping justification
		if self.is_not_null():
			for mol_state in self._excite_string:
				if mol_state==0:  string += ". "
				else:             string += "{} ".format(mol_state)
			string += ">"
		else:
			string += "NULL VECTOR >"
		return string
	def value(self):
		if self.is_not_null():  return (copy.copy(self._excite_string),self._coefficient)
		else:                   raise  Exception("attempted access of occupation information for Null configuration")
	def transition(self, molecule, initial, final):
		""" operates on configuration in place! """
		if self.is_not_null():			# Action onto Null leaves Null.
			if self._excite_string[molecule]==initial:
				self._excite_string[molecule] = final
			else:
				self._make_null()	# If not in appropriate initial state, then it becomes Null.
	def is_not_null(self):
		return (self._excite_string is not None)
	def _make_null(self):
		self._excited_string = None
		self._coefficient     = 0
	def _op_string_action(self,ops):
		new_config = copy.deepcopy(self)
		for op in reversed(ops.operators):  new_config.transition(*op.arguments)
		return new_config
	def __or__(self,other):
		if self._excite_string==other._excite_string:  return self._coefficient * other._coefficient
		else:  return 0




class transition(object):
	def __init__(self, *arguments):
		self.arguments = arguments
	def __call__(self, the_config):
		molecule, initial, final = self.arguments
		new_config = copy.deepcopy(the_config)
		new_config.transition(molecule, initial, final)
		return new_config

class op_string(object):
	def __init__(self,*operators):
		self.operators = operators
	def __or__(self,config_or_state):
		return config_or_state._op_string_action(self)
