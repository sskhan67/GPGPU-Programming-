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

def commute_composite( hamiltonian, excitation_op, tools ):
	com1    = tools.rule_distributor(hamiltonian.right_operator, excitation_op, tools)
	com2    = tools.rule_distributor(hamiltonian.left_operator,  excitation_op, tools)
	num_hamiltonian     = len(com1)
	num_excitation_op   = len(com2)
	resulting_operators = []
	for i in range(num_hamiltonian):
		new_composite = tools.composite_operator( +1, hamiltonian.left_operator,   com1[i]    )
		resulting_operators.append( new_composite )
	for j in range(num_excitation_op):
		new_composite = tools.composite_operator(   +1, com2[j],   hamiltonian.right_operator )
		resulting_operators.append( new_composite )
	return( resulting_operators )

def commute_one_body_one_body(hamiltonian, excitation_op, operator):
	result      = []
	new_offset  = hamiltonian.offset + hamiltonian.offset
	new_coeff   = hamiltonian.coeff * excitation_op.coeff
	if hamiltonian.indices[1] == excitation_op.indices[0]:
		new_hamiltonian     = [+ new_coeff, new_offset, [hamiltonian.indices[0], excitation_op.indices[1] ] ]
		one_bdy_hamiltonian = operator(new_hamiltonian)
		result.append(one_bdy_hamiltonian)
	if hamiltonian.indices[0] == excitation_op.indices[1]:
		new_excitation_op = [ - new_coeff, new_offset, [ excitation_op.indices[0], hamiltonian.indices[1] ] ]
		one_bdy_excitation_op = operator(new_excitation_op)
		result.append(one_bdy_excitation_op)
	return( result )

def commute_two_body_one_body(hamiltonian, excitation_op, operator):
	result    = []
	new_coeff = hamiltonian.coeff * excitation_op.coeff
	if hamiltonian.indices[0] == excitation_op.indices[1]:
		w1          = hamiltonian.offset[0] + excitation_op.offset
		new_offset  = [w1, hamiltonian.offset[1], hamiltonian.offset[2] ]
		new_indices = [excitation_op.indices[0], hamiltonian.indices[1], hamiltonian.indices[2], hamiltonian.indices[3] ]
		new_hamiltonian     = [- new_coeff , new_offset, new_indices ]
		two_bdy_hamiltonian = operator(new_hamiltonian)
		result.append( two_bdy_hamiltonian )
	if hamiltonian.indices[3] == excitation_op.indices[0]:
		w3          = hamiltonian.offset[2] - excitation_op.offset
		new_offset  = [ hamiltonian.offset[0], hamiltonian.offset[1], w3 ]
		new_indices = [ hamiltonian.indices[0], hamiltonian.indices[1], hamiltonian.indices[2], excitation_op.indices[1] ]
		new_excitation_op     = [+1, new_offset, new_indices ]
		two_bdy_excitation_op = operator(new_excitation_op)
		result.append( two_bdy_excitation_op )
	if hamiltonian.indices[1] == excitation_op.indices[1]:
		w2          = hamiltonian.offset[1] + excitation_op.offset
		new_offset  = [ hamiltonian.offset[0],w2 ,hamiltonian.offset[2] ]
		new_indices = [ hamiltonian.indices[0], excitation_op.indices[0], hamiltonian.indices[2], hamiltonian.indices[3] ]
		new_op3     = [-1, new_offset, new_indices ]
		two_bdy_op3 = operator(new_op3)
		result.append( two_bdy_op3 )
	if hamiltonian.indices[2] == excitation_op.indices[0]:
		w1          = hamiltonian.offset[0] + excitation_op.offset
		w2          = hamiltonian.offset[1] + excitation_op.offset
		w4          = hamiltonian.offset[2] + excitation_op.offset
		new_offset  = [ w1, w2, w4]
		new_indices = [ hamiltonian.indices[0], hamiltonian.indices[1], excitation_op.indices[1], hamiltonian.indices[3] ]
		new_op4     = [+1, new_offset, new_indices ]
		two_bdy_op4 = operator(new_op4)
		result.append( two_bdy_op4 )
	return( result )


def commute_one_body_two_body(hamiltonian, excitation_op, operator):
	result = []
	new_coeff = hamiltonian.coeff * excitation_op.coeff
	if hamiltonian.indices[1] == excitation_op.indices[0]:
		w1          = hamiltonian.offset + excitation_op.offset[0]
		new_offset  = [w1, excitation_op.offset[1], excitation_op.offset[2] ]
		new_indices = [hamiltonian.indices[0], excitation_op.indices[1], excitation_op.indices[2], op2.indices[3] ]
		new_hamiltonian     = [+ new_coeff, new_offset, new_indices ]
		two_bdy_hamiltonian = operator(new_hamiltonian)
		result.append(two_bdy_hamiltonian)
	if hamiltonian.indices[0] == excitation_op.indices[3]:
		w2          = excitation_op.offset[2] - hamiltonian.offset
		new_offset  = [ excitation_op.offset[0], excitation_op.offset[1], w2 ]
		new_indices = [ excitation_op.indices[0], excitation_op.indices[1], op2.indices[2], hamiltonian.indices[1] ]
		new_excitation_op     = [- new_coeff, new_offset, new_indices ]
		two_bdy_excitation_op = operator(new_excitation_op)
		result.append(two_bdy_excitation_op)
	if hamiltonian.indices[1] == excitation_op.indices[1]:
		w3          = excitation_op.offset[1] + hamiltonian.offset
		new_offset  = [ excitation_op.offset[0],w3 ,excitation_op.offset[2] ]
		new_indices = [ excitation_op.indices[0], hamiltonian.indices[0], excitation_op.indices[2], op2.indices[3] ]
		new_op3     = [+ new_coeff, new_offset, new_indices ]
		two_bdy_op3 = operator(new_op3)
		result.append(two_bdy_op3)
	if hamiltonian.indices[0] == excitation_op.indices[2]:
		w1          = hamiltonian.offset + excitation_op.offset[0]
		w3          = excitation_op.offset[1] + hamiltonian.offset
		w4          = excitation_op.offset[2] + hamiltonian.offset
		new_offset  = [ w1, w3, w4]
		new_indices = [ excitation_op.indices[0], excitation_op.indices[1], hamiltonian.indices[1], op2.indices[3] ]
		new_op4     = [- new_coeff, new_offset, new_indices ]
		two_bdy_op4 = operator(new_op4)
		result.append(two_bdy_op4)
	return( result )


def commute_two_body_two_body(hamiltonian, excitation_op, tools):
	result    = []
	new_coeff = hamiltonian.coeff * excitation_op.coeff
	if hamiltonian.indices[2] == excitation_op.indices[0]:
		w1           = hamiltonian.offset[0] + excitation_op.offset[0] - excitation_op.offset[2]
		w2           = hamiltonian.offset[1] + excitation_op.offset[0] - excitation_op.offset[2]
		w3           = hamiltonian.offset[2] + excitation_op.offset[0] - excitation_op.offset[2]
		new_offsets  = [w1, w2, w3]
		new_2bd_indx = [hamiltonian.indices[0], hamiltonian.indices[1], excitation_op.indices[3], hamiltonian.indices[3] ]
		new_1bd_indx = [ excitation_op.indices[1], excitation_op.indices[2] ]
		new_2bd_hamiltonian  = [ +1, new_offsets,   new_2bd_indx  ]
		new_1bd_hamiltonian  = [ +1, excitation_op.offset[1], new_1bd_indx ]
		two_bdy_hamiltonian  = tools.operator(new_2bd_hamiltonian)
		one_bdy_hamiltonian  = tools.operator(new_1bd_hamiltonian)
		compost_hamiltonian  = tools.composite_operator( -new_coeff, two_bdy_hamiltonian, one_bdy_hamiltonian ) 
		result.append(compost_hamiltonian) 
	if hamiltonian.indices[3] == excitation_op.indices[0]:
		w1            = excitation_op.offset[0] - hamiltonian.offset[2] + hamiltonian.offset[0]
		w2            = excitation_op.offset[0] - hamiltonian.offset[2] + hamiltonian.offset[1]
		w3            = excitation_op.offset[0] - hamiltonian.offset[2]
		_2bdy_offsets = [w1, w2, w3]
		_2bdy_indx    = [hamiltonian.indices[0], hamiltonian.indices[1], excitation_op.indices[3], hamiltonian.indices[2] ]
		_1bdy_indx    = [excitation_op.indices[1], excitation_op.indices[2] ]
		_2bdy_excitation_op     = [+1, _2bdy_offsets, _2bdy_indx ]
		_1bdy_excitation_op     = [+1, excitation_op.offset[1], _1bdy_indx ]
		two_bdy_excitation_op   = tools.operator(_2bdy_excitation_op)
		one_bdy_excitation_op   = tools.operator(_1bdy_excitation_op)
		compost_excitation_op   = tools.composite_operator( +new_coeff, two_bdy_excitation_op, one_bdy_op2 )
		result.append(compost_excitation_op)
	if hamiltonian.indices[2] == excitation_op.indices[1]:
		w1            = excitation_op.offset[1] - excitation_op.offset[2] + hamiltonian.offset[0]
		w2            = excitation_op.offset[1] - excitation_op.offset[2] + hamiltonian.offset[1]
		w3            = excitation_op.offset[1] - excitation_op.offset[2] + hamiltonian.offset[2]
		_2bdy_offsets = [w1, w2, w3]
		_2bdy_indx    = [hamiltonian.indices[0], hamiltonian.indices[1], excitation_op.indices[3], hamiltonian.indices[3]]
		_1bdy_indx    = [excitation_op.indices[0], excitation_op.indices[2] ]
		_2bdy_op3     = [+1, _2bdy_offsets, _2bdy_indx]
		_1bdy_op3     = [+1, excitation_op.offset[0], _1bdy_indx]
		two_bdy_op3   = tools.operator(_2bdy_op3)
		one_bdy_op3   = tools.operator(_1bdy_op3)
		compost_op3   = tools.composite_operator(+new_coeff, two_bdy_op3, one_bdy_op3)
		result.append(compost_op3)
	if hamiltonian.indices[3] == excitation_op.indices[1]:
		w1            =  excitation_op.offset[1] - hamiltonian.offset[2] - excitation_op.offset[2] + hamiltonian.offset[0]
		w2            =  excitation_op.offset[1] - hamiltonian.offset[2] - excitation_op.offset[2] + hamiltonian.offset[1]
		w3            =  excitation_op.offset[1] - hamiltonian.offset[2] - excitation_op.offset[2]
		_2bdy_offsets = [w1, w2, w3]
		_2bdy_indx    = [hamiltonian.indices[0], hamiltonian.indices[1], excitation_op.indices[3], hamiltonian.indices[2] ]
		_1bdy_indx    = [excitation_op.indices[0], excitation_op.indices[2] ]
		_2bdy_op4     = [+1, _2bdy_offsets, _2bdy_indx]
		_1bdy_op4     = [+1, excitation_op.offset[0], _1bdy_indx]
		two_bdy_op4   = tools.operator(_2bdy_op4)
		one_bdy_op4   = tools.operator(_1bdy_op4) 
		compost_op4   = tools.composite_operator(-new_coeff, two_bdy_op4, one_bdy_op4)
		result.append(compost_op4)
	if hamiltonian.indices[1] == excitation_op.indices[3]:
		w1            =  excitation_op.offset[0] 
		w2            =  excitation_op.offset[1] 
		w3            =  excitation_op.offset[2] - hamiltonian.offset[1] + hamiltonian.offset[2]
		_2bdy_offsets = [w1, w2, w3]
		_2bdy_indx    = [excitation_op.indices[0], excitation_op.indices[1], op2.indices[2], hamiltonian.indices[3] ]
		_1bdy_indx    = [hamiltonian.indices[0], hamiltonian.indices[2] ]
		_2bdy_op5     = [+1, _2bdy_offsets, _2bdy_indx]
		_1bdy_op5     = [+1, hamiltonian.offset[0], _1bdy_indx]
		two_bdy_op5   = tools.operator(_2bdy_op5)
		one_bdy_op5   = tools.operator(_1bdy_op5)
		compost_op5   = tools.composite_operator(-new_coeff, two_bdy_op5, one_bdy_op5)
		result.append(compost_op4)
	if hamiltonian.indices[0] == excitation_op.indices[3]:
		w1            =  excitation_op.offset[0]
		w2            =  excitation_op.offset[1]
		w3            =  excitation_op.offset[2] - hamiltonian.offset[0] + hamiltonian.offset[2]
		_2bdy_offsets = [w1, w2, w3]
		_2bdy_indx    = [excitation_op.indices[0], excitation_op.indices[1], op2.indices[2], hamiltonian.indices[3] ]
		_1bdy_indx    = [hamiltonian.indices[1], hamiltonian.indices[2] ]
		_2bdy_op6     = [+1, _2bdy_offsets, _2bdy_indx]
		_1bdy_op6     = [+1, hamiltonian.offset[1], _1bdy_indx]
		two_bdy_op6   = tools.operator(_2bdy_op6)
		one_bdy_op6   = tools.operator(_1bdy_op6)
		compost_op6   = tools.composite_operator(+new_coeff, two_bdy_op6, one_bdy_op6)
		result.append(compost_op6)
	if hamiltonian.indices[0] == excitation_op.indices[2]:
		w1            =  hamiltonian.offset[0] - hamiltonian.offset[2] + excitation_op.offset[0]
		w2            =  hamiltonian.offset[0] - hamiltonian.offset[2] + excitation_op.offset[1]
		w3            =  hamiltonian.offset[0] - hamiltonian.offset[2] + excitation_op.offset[2]
		_2bdy_offsets = [w1, w2, w3]
		_2bdy_indx    = [excitation_op.indices[0], excitation_op.indices[1], hamiltonian.indices[3], op2.indices[3] ]
		_1bdy_indx    = [hamiltonian.indices[1], hamiltonian.indices[2] ]
		_2bdy_op7     = [+1, _2bdy_offsets, _2bdy_indx]
		_1bdy_op7     = [+1, hamiltonian.offset[1], _1bdy_indx]
		two_bdy_op7   = tools.operator(_2bdy_op7)
		one_bdy_op7   = tools.operator(_1bdy_op7)
		compost_op7   = tools.composite_operator(+new_coeff, two_bdy_op7, one_bdy_op7)
		result.append(compost_op7)
	if hamiltonian.indices[1] == excitation_op.indices[2]:
		w1            =  hamiltonian.offset[1] - hamiltonian.offset[2] + excitation_op.offset[0]
		w2            =  hamiltonian.offset[1] - hamiltonian.offset[2] + excitation_op.offset[1]
		w3            =  hamiltonian.offset[1] - hamiltonian.offset[2] + excitation_op.offset[2]
		_2bdy_offsets = [w1, w2, w3]
		_2bdy_indx    = [excitation_op.indices[0], excitation_op.indices[1], hamiltonian.indices[3], op2.indices[3] ]
		_1bdy_indx    = [hamiltonian.indices[0], hamiltonian.indices[2] ]
		_2bdy_op8     = [+1, _2bdy_offsets, _2bdy_indx]
		_1bdy_op8     = [+1, hamiltonian.offset[0], _1bdy_indx]
		two_bdy_op8   = tools.operator(_2bdy_op8)
		one_bdy_op8   = tools.operator(_1bdy_op8)
		compost_op8   = tools.composite_operator(-new_coeff, two_bdy_op8, one_bdy_op8)
		result.append(compost_op8)
	return( result )
