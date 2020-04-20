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
# This file tests the commutator by creating operarators and commuting them.

from   classes          import operator, composite_operator, commutator, tools_encapsulation
from   function         import rule_distributor
from   parameters       import N, c1, e, f, d, p, r, q, s
from   parameters       import c2, k, l, h, a, b, i, j
from   primitive_rules  import commute_one_body_one_body, commute_two_body_two_body, commute_one_body_two_body, commute_two_body_one_body


tools = tools_encapsulation( operator, composite_operator, rule_distributor, commute_one_body_one_body, 
       commute_two_body_two_body, commute_one_body_two_body, commute_two_body_one_body  )


# This block initiatate the hamiltonian two and one body operators, 
# along with the excitation operators 
 
two_body_operator_parameters  = [+c1, [e, f, d], [p, r, q, s] ]
double_excitation_parameters  = [ c2, [k, l, h], [a, b, i, j] ]
one_body_operator_parameters  = [+c1, e, [p, q] ]
single_exitation_parameters   = [ c2, k, [a, i] ] 
two_body_operator             = operator( two_body_operator_parameters )
double_excitation_operator    = operator( double_excitation_parameters )



# this block commutes the relevant operators

commute_operators   = commutator( tools )
resulting_operators = commute_operators(two_body_operator, double_excitation_operator)

print( resulting_operators[0].full_op )




'''
A = 'A'
B = 'B'
I = 'I'
J = 'J'

T2  = [ c2, [k, l, h], [A, B, I, J] ]
op3 = operator(T2)

com3 = com1(com2[0], op3)
print(com3[0].full_op)
'''
