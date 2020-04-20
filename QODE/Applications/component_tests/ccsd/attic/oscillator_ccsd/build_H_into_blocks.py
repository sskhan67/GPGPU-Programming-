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
#
# ANALYTICAL SOLUTION = 2.398435656980412
# CISD Results:
#              5 x 5 => 2.3984416499
#              7 x 7 => 2.39843571482
#

import sys
import time
import copy
import numpy as np
# from qode.util import parallel, output, textlog, indent

# Config and state infrastructure, SCF and CCSD modules
# from qode.fermion_field import state
# from qode.SCF_diag_full import hf as scf_routine
# from qode.many_body     import CCSD
from Applications.component_tests.oscillator_ccsd.old_qode.generate_ho_config import diag_rec_h0_mat, diag_rec_v_matrix, generate_main, get_lowest_analytical, get_rec_obj

# Oscillator System Imports
from qode.coupled_oscillators            import oscillator_system as osys
from qode.coupled_oscillators            import hamiltonian_class
from qode.coupled_oscillators.analytical import analytical_solution
from qode.many_body.fast_operator.fast_operator import cc_operator
from qode.many_body.fast_operator.generalutils  import print_array 

from copy import deepcopy


# def init_H_blocks_zero(cc_op_obj, No, Nv):
# 	# cc_op_obj.K      = np.zeros(1)        # This is zero
# 	# cc_op_obj.Ex     = np.zeros((No*Nv))  # This is zero
# 	cc_op_obj.Fo     = np.zeros((No*No))
# 	cc_op_obj.Fv     = np.zeros((Nv*Nv))
# 	# cc_op_obj.Dx     = np.zeros((Nv*No))  # This is zero
# 	cc_op_obj.ExEx   = np.zeros(((No*Nv)**2))
# 	cc_op_obj.ExFo   = np.zeros((Nv*No**3))
# 	cc_op_obj.ExFv   = np.zeros((No*Nv**3))
# 	cc_op_obj.ExDx   = np.zeros(((No*Nv)**2))
# 	cc_op_obj.FoFo   = np.zeros((No**4))
# 	cc_op_obj.FvFv   = np.zeros((Nv**4))
# 	cc_op_obj.FoDx   = np.zeros((Nv*No**3))
# 	cc_op_obj.FvDx   = np.zeros((No*Nv**3))
# 	cc_op_obj.DxDx   = np.zeros(((No*Nv)**2))	
# 	return cc_op_obj

# def abs_to_relative_pos(this_index, occ_orbs, vrt_orbs):
# 	No = len(occ_orbs)
# 	Nv = len(vrt_orbs)
# 	# print(occ_orbs,vrt_orbs,this_index)
# 	if this_index in vrt_orbs:
# 		return vrt_orbs.index(this_index), Nv, True
# 	elif this_index in occ_orbs:
# 		return occ_orbs.index(this_index), No, False
# 	else:
# 		raise AssertionError


# def fill_in_h_mat(H, diag_index, value, occ_orbs, vrt_orbs):
# 	# h_mat is STRICTLY DIAGONAL, SAVE SOME TIME...
# 	this_index, multiplier, ex_flag = abs_to_relative_pos(diag_index, occ_orbs, vrt_orbs)
# 	# print(this_index,multiplier,ex_flag)
# 	if ex_flag:  # True => Virtual
# 		H.Fv[ this_index * (1 + multiplier) ] += value
# 	else:
# 		H.Fo[ this_index * (1 + multiplier) ] += value

	

# dict_H = { 'ExEx':[False,True,False,True], 'ExFo':[False,True,False,False] , 'ExFv':[False,True,True,True], \
#            'ExDx':[False,True,True,False], 'FoFo':[False,False,False,False], 'FvFv':[True,True,True,True],  \
#            'FoDx':[False,False,True,False], 'FvDx':[True,True,True,False]   , 'DxDx':[True,False,True,False] } 

# def fill_in_V_mat(H, V_mat_obj, mol_pair, indices, rec_num_states):
# 	inputs = mol_pair + indices
# 	value  = V_mat_obj(*inputs)
# 	block_indices = []
# 	for item in indices:
# 		if item == 0:
# 			block_indices += [False]
# 		else:
# 			block_indices += [True]
# 	# print(mol_pair, indices, value, block_indices)
# 	block_key = False
# 	for key, block_flag in dict_H.items():
# 		if block_flag == block_indices:
# 			block_key = key
# 	# print(mol_pair, indices, value, block_indices, block_key)
# 	if block_key != False:
# 		# print(mol_pair, indices, value, block_indices, block_key)
# 		real_indices = deepcopy(indices)
# 		x, y = mol_pair
# 		for i in range(x):
# 			real_indices[0] += rec_num_states[i]
# 			real_indices[2] += rec_num_states[i]
# 			if block_indices[0] == True:
# 				real_indices[0] += 1
# 			if block_indices[2] == True:
# 				real_indices[2] += 1
# 		for i in range(y):
# 			real_indices[1] += rec_num_states[i]
# 			real_indices[3] += rec_num_states[i]
# 			if block_indices[1] == True:
# 				real_indices[1] += 1
# 			if block_indices[3] == True:
# 				real_indices[3] += 1
# 		print(mol_pair, indices, real_indices, value, block_indices, block_key)
# 		H.__dict__[block_key] += value

# def build_H(rec_num_states, h_mat_obj, V_mat_obj):
# 	No = len(rec_num_states)
# 	Nv = 0
# 	for item in rec_num_states:
# 		Nv += (item - 1)
# 	#
# 	occ_orbs = [0]
# 	vrt_orbs = [ occ_orbs[-1]+1+i for i in range(rec_num_states[0]-1) ]
# 	for n in rec_num_states[1:]:
# 		occ_orbs += [ vrt_orbs[-1]+1 ]
# 		vrt_orbs += [ occ_orbs[-1]+1+i for i in range(n-1) ]
# 	print("OCC =", No, occ_orbs)
# 	print("VRT =", Nv, vrt_orbs)
# 	H = init_H_blocks_zero(cc_operator(), No, Nv)
# 	# print_array(H)
# 	for i in range(h_mat_obj.get_dimension()):
# 		fill_in_h_mat(H, i, h_mat_obj(i,i), occ_orbs, vrt_orbs )

# 	num_mol = len(rec_num_states)
# 	print("NUM MOL =",num_mol)
# 	for x in range(num_mol):
# 		for y in range(x+1, num_mol):
# 			# V Pairs now!
# 			for i in range(rec_num_states[x]):
# 				for j in range(rec_num_states[y]):
# 					for k in range(i+1, rec_num_states[x]):
# 						for l in range(j+1, rec_num_states[y]):
# 							# V_{i,j,k,l} is being filled into double-operator blocks
# 							fill_in_V_mat(H, V_mat_obj, [x,y], [i,j,k,l], rec_num_states)


# 	print_array(H)





if __name__ == "__main__":
	# out = output(log=textlog(echo=True))
	# resources = parallel.resources(int(sys.argv[1]))
	# 
	rec_num_states = [11, 11]
	
	ho1 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho2 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	

	ho3 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho4 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	
	
	coupling_mat = [[0.0, 0.3],\
	                [0.3, 0.0]]

	
	mol1 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho1,ho2] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	
	mol2 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho3,ho4] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	
	
	ext_k_mat = [[0.0,0.1],[0.1,0.0]]
	
	mol_obj = osys.OscillatorSystem( osys.group_of_oscillators( recursive_list =  [mol1, mol2],
                                                                 coupling     =  ext_k_mat ) )
	
	rec_obj = get_rec_obj( mol_obj, rec_num_states )
	
	
	eig1, eig2 = rec_obj.get_sorted_eigen_energy()
	v_mat_list = rec_obj.get_pair_dipole_mat()
	

	# Building h_mat and V_mat objects
	h_mat = diag_rec_h0_mat([eig1, eig2])
	V_mat = diag_rec_v_matrix(rec_num_states, v_mat_list)


	print("H original dimension =", h_mat.get_dimension())
	build_H(rec_num_states, h_mat, V_mat)




