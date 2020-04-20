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
import Applications.general_ci.config as config
import Applications.general_ci.state  as state
from Applications.osc_explicit_state.generate_ho_config import generate_singles


def T1(state_obj, t_amp_obj, rec_num_states):
	# Take any state_obj but make a new_state
	#
	# T1 = \sum_{a i} t_ia a+ i
	# 
	#   t_ia are the coeffs in state_obj but only the SINGLE excitation amplitudes.
	#
	new_state = state.state(state_obj.get_configs())
	t_amps    = t_amp_obj.get_coeffs()
	idx_shift = 0
	num_mol   = len(rec_num_states)
	t_idx     = 1
	for i in range(num_mol):
		for j in range(1, rec_num_states[i]):
			# print('mol =',i,'ref =>',j)
			#
			#
			op_list     = [[config.create, idx_shift + j ], [config.annihilate, idx_shift]]
			temp_state  = state.operate(op_list, state_obj)
			new_state.add_to( t_amps[ t_idx ] , temp_state )
			t_idx += 1
		#		
		idx_shift += rec_num_states[i]

	return new_state




def T2(state_obj, t_amp_obj, rec_num_states):
	# Take any state_obj but make a new_state
	#
	# T2 = \sum_{a<b, i<j} t_{ijab} a+ i b+ j 
	#
	#   t_{ijab} are the coeffs in state_obj but only the DOUBLE excitation amplitudes.
	#
	new_state = state.state(state_obj.get_configs())
	t_amps    = t_amp_obj.get_coeffs()	
	#
	num_mol   = len(rec_num_states)
	#
	#
	t_idx = 1 + len(generate_singles(rec_num_states))
	first_idx_shift = 0
	for i in range(num_mol):   # Mol 1
		#
		#
		second_idx_shift = 0
		for num in rec_num_states[:i+1]:
			second_idx_shift += num
		#
		#
		for j in range(i+1, num_mol): # Mol 2
			#
			#
			for a in range(1, rec_num_states[i]):
				#
				#
				for b in range(1, rec_num_states[j]):
					#
					#
					op_list    = [[config.create, first_idx_shift + a], [config.annihilate, first_idx_shift], [config.create, second_idx_shift + b], [config.annihilate, second_idx_shift] ]
					temp_state = state.operate(op_list, state_obj)
					new_state.add_to( t_amps[ t_idx ], temp_state )
					t_idx += 1
			#
			#
			second_idx_shift += rec_num_states[j]
		first_idx_shift += rec_num_states[i]
	#
	#
	return new_state




def T(state_obj, t_amp_obj, h_mat, V_mat, rec_num_states):
	return state.add_state( T1(state_obj, t_amp_obj, rec_num_states),\
							T2(state_obj, t_amp_obj, rec_num_states) )







