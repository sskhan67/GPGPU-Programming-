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
import numpy  as np
import ctypes as ct
import copy

class Storage(ct.Structure):
	# each struct stores one operator sequence storgae
	_fields_ = [("name", ct.c_char * 5),
                ("data", ct.POINTER(ct.POINTER(ct.POINTER(ct.c_double))))]

class ContractedHamiltonian(ct.Structure):      
	# each struct stores part of the Hamiltonian for fixed number of electrons
	_fields_ = [("num_elec", ct.c_ssize_t),
				("pairs",    ct.POINTER(ct.c_ssize_t)),
				("n_pairs",  ct.c_ssize_t),
				("data",     ct.POINTER(ct.POINTER(ct.c_double)))]

def parse_num_states_used(dense_mat, n_ch_st):
	num_states = []  # num states used
	for i in range(n_ch_st):
		num_states.append( len(dense_mat['ca'][i][i]) )
	return num_states

def parse_num_spin_orbs(dense_mat, n_ch_st):
	return dense_mat['ca'][0][0][0][0].shape[0]

def operator_to_c_struct(op_storage, n_ch_st, num_states):
	#
	dim = n_ch_st * n_ch_st
	p_op_storage = (ct.POINTER(ct.POINTER(ct.c_double)) * dim )()
	#
	for i in range(n_ch_st):
		for j in range(n_ch_st):
			if op_storage[i][j]:
				sub_dim = num_states[i] * num_states[j]
				p_op_storage[i*n_ch_st+j] = (ct.POINTER(ct.c_double) * sub_dim)()
				for k in range(num_states[i]):
					for l in range(num_states[j]):
						p_op_storage[i*n_ch_st+j][k*num_states[j]+l] = op_storage[i][j][k][l].ctypes.data_as( ct.POINTER(ct.c_double) )
			else:
				p_op_storage[i*n_ch_st+j] = None
	return p_op_storage

def dense_mat_to_c_struct(dense_mat):
	#
	op_types   = ['a', 'c', 'aa', 'cc', 'ca', 'caa', 'cca', 'ccaa']
	n_ch_st    = len( dense_mat['a'] )
	num_states = parse_num_states_used(dense_mat, n_ch_st)
	C_storage  = (Storage * 8)() 
	count = 0
	for key in op_types:
		C_storage[count].name = key.encode('ascii')
		C_storage[count].data = operator_to_c_struct(dense_mat[key], n_ch_st, num_states)
		count += 1
	return C_storage, n_ch_st, num_states

def two_body_num_elec_pairs(num_elec1, num_elec2): # num_elec is naturally sorted
	lo1 = min(num_elec1)
	hi1 = max(num_elec1)
	lo2 = min(num_elec2)
	hi2 = max(num_elec2)
	num_elec_two_body = []
	possible_pairs = []
	for n_elec in range(lo1+lo2, hi1+hi2+1):
		tmp = []
		for n1 in num_elec1: # This is very generic, it fits discrete numbers of electrons
			n2 = n_elec - n1
			if n2 in num_elec2:
				tmp.append( (n1,n2) )
		if tmp:
			num_elec_two_body.append(n_elec)
			possible_pairs.append( list(tmp) )
	return possible_pairs, num_elec_two_body

def create_two_body_Ham_storage(num_states1, num_states2, num_elec1, num_elec2, num_elec_pairs):
	H = []
	index_list = []
	for item in (num_elec_pairs):
		tmp = []
		tmp_index = []
		for bra in item:
			i, j  = bra
			i_idx,     j_idx     = num_elec1.index(i), num_elec2.index(j)
			n_state_i, n_state_j = num_states1[i_idx], num_states2[j_idx]
			tmp_index.append( [i_idx,j_idx] )
			for ket in item:
				k, l = ket
				k_idx,     l_idx     = num_elec1.index(k), num_elec2.index(l)
				n_state_k, n_state_l = num_states1[k_idx], num_states2[l_idx]
				tmp.append( np.zeros((n_state_i*n_state_j, n_state_k*n_state_l)) )
		H.append( copy.deepcopy(tmp) )
		index_list.append( list(tmp_index) )
	return H, index_list

def fixed_n_elec_Ham_to_c_struct(fixed_e_H):
	n = len(fixed_e_H)
	C_pp_double = (ct.POINTER(ct.c_double) * n)()
	for i in range(n):
		C_pp_double[i] = fixed_e_H[i].ctypes.data_as( ct.POINTER(ct.c_double) )
	return C_pp_double

def Hamiltonian_to_c_struct(Ham, num_elec_two_body, index_list):
	n = len(Ham)
	C_contHam = (ContractedHamiltonian * n)()
	for i in range(n):
		C_contHam[i].num_elec = num_elec_two_body[i]
		C_contHam[i].n_pairs  = len(index_list[i])
		dim = 2 * len(index_list[i])
		C_contHam[i].pairs    = (ct.c_ssize_t * dim)()
		C_contHam[i].data     = fixed_n_elec_Ham_to_c_struct(Ham[i])
		count = 0
		for j in range(len(index_list[i])):
			for k in range(2):
				C_contHam[i].pairs[count] = index_list[i][j][k]
				count += 1
	return C_contHam


def main(num_elec, dense_mat, hmat, Vmat):
	C_store, n_ch_st, num_states = dense_mat_to_c_struct(dense_mat)
	
	n_spin = parse_num_spin_orbs(dense_mat, n_ch_st)
	
	C_num_states = (ct.c_ssize_t * n_ch_st)()
	for i in range(n_ch_st): C_num_states[i] = num_states[i]

	C_num_elec = (ct.c_ssize_t * n_ch_st)()
	for i in range(n_ch_st): C_num_elec[i]   = num_elec[i]

	C_hmat = hmat.ctypes.data_as( ct.POINTER(ct.c_double) )
	C_Vmat = Vmat.ctypes.data_as( ct.POINTER(ct.c_double) )

	num_elec_pairs, num_elec_two_body = two_body_num_elec_pairs(num_elec, num_elec)

	H, index_list = create_two_body_Ham_storage(num_states, num_states, num_elec, num_elec, num_elec_pairs)
	C_contH       = Hamiltonian_to_c_struct(H, num_elec_two_body, index_list)
	num_H_mat     = len(H)

	lib = ct.cdll.LoadLibrary('./numerical_engines/libbiorthohamcontract.so')

	lib.contract_biorthogonal_ham( ct.c_ssize_t(n_ch_st), ct.c_ssize_t(n_spin), C_num_states, C_num_elec,
                                   ct.c_ssize_t(n_ch_st), ct.c_ssize_t(n_spin), C_num_states, C_num_elec,
                                   C_store, C_store,
                                   C_hmat, C_Vmat, 
                                   C_contH, ct.c_ssize_t(num_H_mat), ct.c_ssize_t(1) )

	return H


if __name__ == '__main__':
	import sys
	import pickle

	num_elec  = [3, 4, 5]
	dense_mat = pickle.load( open('dense_mat.p','rb') )
	hmat = np.load('biorthogonal_core_matrix_dimer.npy')
	Vmat = np.load('biorthogonal_twobody_matrix_dimer.npy')


	H_mats = main(num_elec, dense_mat, hmat, Vmat)

	# for item in H_mats:
	# 	print('-'*80)
	# 	for mat in item:
	# 		print("max = %.3e min = %.3e" %(mat.max(), mat.min()))
	# print('-'*80)


	pickle.dump(H_mats, open('H.p','wb'))

	# for item in num_elec_pairs:
	# 	print(item)
	# print(num_elec_two_body)
	# for item in index_list:
	# 	print(item)

	# for item in H:
	# 	print(len(item))
	# 	for mat in item:
	# 		print(mat.shape)

