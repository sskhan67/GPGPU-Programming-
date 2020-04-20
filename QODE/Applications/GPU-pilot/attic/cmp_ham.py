#    (C) Copyright 2018 Yuhong Liu
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
import os
import sys
import pickle
import numpy as np
from biorthogonal_Hamiltonian import two_body_num_elec_pairs 
# def two_body_num_elec_pairs(num_elec1, num_elec2)



def tony_to_sparse(V, num_elec, num_states_used):
	tony_order = []
	for item1 in num_elec:
		for item2 in num_elec:
			tony_order.append([item1,item2])
	V_dict = { }
	dim_a = 0
	for i,j in tony_order:
		idx_i, idx_j = num_elec.index(i), num_elec.index(j)
		dim1, dim2 = num_states_used[idx_i], num_states_used[idx_j]
		tmp_dim_a = (dim_a, dim_a + dim1 * dim2)
		dim_b = 0
		for k,l in tony_order:
			idx_k, idx_l = num_elec.index(k), num_elec.index(l)
			dim3, dim4 = num_states_used[idx_k], num_states_used[idx_l]
			tmp_dim_b = (dim_b, dim_b + dim3 * dim4)
			key = i,j,k,l
			V_dict[key] = V[tmp_dim_a[0]:tmp_dim_a[1], tmp_dim_b[0]:tmp_dim_b[1]].copy()
			# print(key, V_dict[key].shape)
			dim_b = tmp_dim_b[1]
		dim_a = tmp_dim_a[1]
	return V_dict


def my_H_to_dict(hmats, my_num_elec):
	pairs, num_elec_two_body = two_body_num_elec_pairs(my_num_elec, my_num_elec)

	keys = []
	for state in pairs:
		tmp = []
		for item1 in state:
			for item2 in state:
				tmp.append( (*item1, *item2) )
		keys.append(list(tmp))
	
	flat_keys = []
	for item in keys:
		flat_keys.extend(item)

	flat_hmat = []
	for item in hmats:
		flat_hmat.extend(item)

	my_V_dict = { }
	for i in range(len(flat_hmat)):
		my_V_dict[flat_keys[i]] = flat_hmat[i]

	return my_V_dict

if __name__ == '__main__':


	V = np.load('V.npy')
	num_elec        = [4,3,5]
	num_states_used = [11, 4, 8]
	Tony_V_dict = tony_to_sparse(V, num_elec, num_states_used)
	print(Tony_V_dict.keys())

	print('-'*80)

	fname = 'H.p'
	hmats = pickle.load( open(fname,'rb') )
	my_num_elec = [3,4,5]
	my_V_dict = my_H_to_dict(hmats, my_num_elec)
	print(my_V_dict.keys())






	for key in my_V_dict.keys():
		my_mat   = my_V_dict[key]
		tony_mat = Tony_V_dict[key]
		if my_mat.shape == tony_mat.shape:
			diff = np.abs( np.abs(my_mat) - np.abs(tony_mat) )
			print(key, diff.max())
		else:
			raise RuntimeError













