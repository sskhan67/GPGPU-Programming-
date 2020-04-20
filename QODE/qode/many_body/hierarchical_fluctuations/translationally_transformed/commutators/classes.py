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


class operator:

	# This class generates one and two body operators and stores its characteristics

	def __init__(self, parameters):

		# This class takes as an argument a list of indices that define the operator
	
		if len(parameters[2]) == 2:
			self.type    = "one body"
			self.coeff   = parameters[0]
			self.offset  = parameters[1]
			self.indices = parameters[2]
			self.full_op = [ self.coeff, self.offset, self.indices ]
		if len(parameters[2]) == 4:
			self.type    = "two body"
			self.coeff   = parameters[0]
			self.offset  = parameters[1]
			self.indices = parameters[2]
			self.full_op = [self.coeff, self.offset, self.indices ]


class composite_operator:

	# This class generates a composite operator that is constituted by two operators
	# A composite operator will often be a composite of composite operators
 
        def __init__(self, coeff, left_operator, right_operator):
                self.coeff             = coeff
                self.type              = "composite_operator"
                self.left_operator     = left_operator
                self.right_operator    = right_operator
                self.full_op           = [self.coeff, self.left_operator.full_op, self.right_operator.full_op]


class commutator:

	# This class takes as arguments the classes to generate its own operators
	# and the function that distributes the commuted operators to their corresponding rules

	def __init__(self, tools):
		self.tools		= tools
		self.rule_distributor	= tools.rule_distributor
	def __call__(self, Hamiltonian_operator, excitation_operator):
		resulting_operators = self.rule_distributor(Hamiltonian_operator, excitation_operator, self.tools)
		return(resulting_operators)

class tools_encapsulation:
	def __init__(self, operator, composite_operator, rule_distributor, commute_one_body_one_body, commute_two_body_two_body, commute_one_body_two_body, commute_two_body_one_body ):
		self.operator                        = operator
		self.composite_operator              = composite_operator
		self.rule_distributor                = rule_distributor
		self.commute_one_body_one_body       = commute_one_body_one_body
		self.commute_two_body_two_body       = commute_two_body_two_body
		self.commute_one_body_two_body       = commute_one_body_two_body
		self.commute_two_body_one_body       = commute_two_body_one_body 
