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
# This Module Loads Molecular Eigenvalues and Paired Coupling Matrix Values into corresponding blocks 
# for Harmonic Oscillator Calculations

# from Applications.reorder_play.test import load_H
# from Applications.reorder_play.maps import map_main
from attic.oscillator_ccsd.old_qode.generate_ho_config import  diag_rec_v_matrix
from qode.many_body.local_operator.local_operator import cc_operator, init_K, init_blocks
import numpy as np

def v_list_find_idx(i, j, j_max):
	if i >= j:
		raise AssertionError
	idx = 0
	for x in range(i):
		idx += j_max-(x+1)-1
	idx += j - 1
	return idx

block_dict ={	'Ex': ['OCC','VRT'],
				'Dx': ['VRT','OCC'],
				'Fo': ['OCC','OCC'],
				'Fv': ['VRT','VRT'],
				'ExEx': ['OCC','VRT','OCC','VRT'],
				'ExFo': ['OCC','VRT','OCC','OCC'],
				'ExFv': ['OCC','VRT','VRT','VRT'],
				'ExDx': ['OCC','VRT','VRT','OCC'],
				'FoFo': ['OCC','OCC','OCC','OCC'],
				'FvFv': ['VRT','VRT','VRT','VRT'],
				'FoDx': ['OCC','OCC','VRT','OCC'],
				'FvDx': ['VRT','VRT','VRT','OCC'],
				'FoFv': ['OCC','OCC','VRT','VRT'],
				'DxDx': ['VRT','OCC','VRT','OCC'],
				'ExExDx': ['OCC','VRT','OCC','VRT','VRT','OCC'],
				'ExFoFo': ['OCC','VRT','OCC','OCC','OCC','OCC'],
				'ExFvFv': ['OCC','VRT','VRT','VRT','VRT','VRT'],
				'ExFoDx': ['OCC','VRT','OCC','OCC','VRT','OCC'],
				'ExFvDx': ['OCC','VRT','VRT','VRT','VRT','OCC'],
				'ExFoFv': ['OCC','VRT','OCC','OCC','VRT','VRT'] }

# def two_idx(i,j):
# 	abs_idx = 0
# 	for x in range(i):
# 		abs_idx += x + 1
# 	abs_idx += j
# 	return abs_idx


# def three_idx(i,j,k):
# 	abs_idx = 0
# 	for x in range(i):
# 		for y in range(x+1):
# 			abs_idx += y + 1
# 	for y in range(j):
# 		abs_idx += y + 1
# 	abs_idx += k
# 	return abs_idx

'''
n_mol = 4
for i in range(n_mol):
	for j in range(i+1, n_mol):
		print("Mol", i,j, "=>", v_list_find_idx(i,j,n_mol))

for m in range(n_mol):
	for n in range(m+1):
		print("Mol", m,n,'->',two_idx(m,n))

for m in range(n_mol):
	for n in range(m+1):
		for k in range(n+1):
			print("Mol",m,n,k,'-->',three_idx(m,n,k))
'''

'''
		self.K      = None
		self.Ex     = None
		self.Fo     = None
		self.Fv     = None
		self.Dx     = None
		self.ExEx   = None   # diag = zero
		self.ExFo   = None   # diag = zero
		self.ExFv   = None   # diag = zero
		self.ExDx   = None
		self.FoFo   = None
		self.FvFv   = None
		self.FoFv   = None   # diag = zero
		# self.FvFo   = None   # diag = zero --> removed in the latest version (Aug 29, 2016)
		self.FoDx   = None   # diag = zero
		self.FvDx   = None   # diag = zero
		self.DxDx   = None
		self.ExFoFo = None
		self.ExFvFv = None
		self.ExFoFv = None
		# self.ExFvFo = None # --> removed in the latest version (Aug 29, 2016)
		self.ExFoDx = None
		self.ExFvDx = None
		self.ExExDx = None '''


def init_H_op(cc_obj, rec):
	cc_obj.K,  cc_obj.K_starters, cc_obj.K_int_dim, cc_obj.K_dim   = init_K()
	cc_obj.Ex, cc_obj.Ex_starters, cc_obj.Ex_int_dim, cc_obj.Ex_dim   = init_blocks(rec, block_dict['Ex'])
	cc_obj.Fo, cc_obj.Fo_starters, cc_obj.Fo_int_dim, cc_obj.Fo_dim   = init_blocks(rec, block_dict['Fo'])
	cc_obj.Fv, cc_obj.Fv_starters, cc_obj.Fv_int_dim, cc_obj.Fv_dim   = init_blocks(rec, block_dict['Fv'])
	cc_obj.Dx, cc_obj.Dx_starters, cc_obj.Dx_int_dim, cc_obj.Dx_dim   = init_blocks(rec, block_dict['Dx'])
	cc_obj.ExEx, cc_obj.ExEx_starters, cc_obj.ExEx_int_dim, cc_obj.ExEx_dim = init_blocks(rec, block_dict['ExEx'])
	cc_obj.ExFo, cc_obj.ExFo_starters, cc_obj.ExFo_int_dim, cc_obj.ExFo_dim = init_blocks(rec, block_dict['ExFo'])
	cc_obj.ExFv, cc_obj.ExFv_starters, cc_obj.ExFv_int_dim, cc_obj.ExFv_dim = init_blocks(rec, block_dict['ExFv'])
	cc_obj.ExDx, cc_obj.ExDx_starters, cc_obj.ExDx_int_dim, cc_obj.ExDx_dim = init_blocks(rec, block_dict['ExDx'])
	cc_obj.FoFo, cc_obj.FoFo_starters, cc_obj.FoFo_int_dim, cc_obj.FoFo_dim = init_blocks(rec, block_dict['FoFo'])
	cc_obj.FvFv, cc_obj.FvFv_starters, cc_obj.FvFv_int_dim, cc_obj.FvFv_dim = init_blocks(rec, block_dict['FvFv'])
	cc_obj.FoDx, cc_obj.FoDx_starters, cc_obj.FoDx_int_dim, cc_obj.FoDx_dim = init_blocks(rec, block_dict['FoDx'])
	cc_obj.FvDx, cc_obj.FvDx_starters, cc_obj.FvDx_int_dim, cc_obj.FvDx_dim = init_blocks(rec, block_dict['FvDx'])
	cc_obj.FoFv, cc_obj.FoFv_starters, cc_obj.FoFv_int_dim, cc_obj.FoFv_dim = init_blocks(rec, block_dict['FoFv'])
	cc_obj.DxDx, cc_obj.DxDx_starters, cc_obj.DxDx_int_dim, cc_obj.DxDx_dim = init_blocks(rec, block_dict['DxDx'])
	return cc_obj

def load_K(H_obj, all_mol_eigenvalues):
	for i in range(len(all_mol_eigenvalues)):
		H_obj.K[0] += all_mol_eigenvalues[i][0]
	return H_obj 


def load_eigen_values_into_H(H_obj, rec, all_mol_eigenvalues):
	num_mol = len(rec)
	Nv = [ num-1 for num in rec ]
	#
	# Fo Parts:  <0|h|0>
	#
	for i in range(num_mol):
		H_obj.Fo[i] = -1.0*all_mol_eigenvalues[i][0]
	#
	# Fv Parts: <a|h|a>
	#
	for n in range(num_mol):
		for a in range(Nv[n]):
			H_obj.Fv[ H_obj.Fv_starters[n] + a * (Nv[n]+1) ] += all_mol_eigenvalues[n][a+1] # This is CORRECT!!
	# Ex Parts: All Zero
	# Dx Parts: All Zero
	return H_obj

def load_v_matrice_into_H(H_obj, rec, v_mat_list):
	num_mol = len(rec)
	Nv = [ num-1 for num in rec ]
	V_mat = diag_rec_v_matrix(rec, v_mat_list)
	#
	#
	# ExEx[i,a,j,b]  m<n: <ij|V|ab>; m>n: <ji|V|ba>
	for m in range(num_mol):
		for n in range(num_mol): 
			if m > n: 
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						# H_obj.ExEx[ H_obj.ExEx_starters[m*num_mol+n] + a*Nv[n] + b ] +=  v_mat_list[ v_list_find_idx(n,m,num_mol) ][0][(b+1)*rec[m] + a+1]
						H_obj.ExEx[ H_obj.ExEx_starters[m*num_mol+n] + a*Nv[n] + b ] +=  V_mat( n, m, 0, 0, b+1, a+1 )


	# ExFo: ExFo[i,a,j,k]  m<n: <ij|V|ak>; m>n: <ji|V|ka>
	for m in range(num_mol):
		for n in range(num_mol): 
			if m < n:
				for a in range(Nv[m]):
					# H_obj.ExFo[ H_obj.ExFo_starters[m*num_mol+n] + a ] +=  -v_mat_list[ v_list_find_idx(m,n,num_mol) ][0][a+1]
					H_obj.ExFo[ H_obj.ExFo_starters[m*num_mol+n] + a ] +=  -V_mat( m, n, 0, 0, a+1, 0 )
			elif m > n:
				for a in range(Nv[m]):
					# H_obj.ExFo[ H_obj.ExFo_starters[m*num_mol+n] + a ] +=  -v_mat_list[ v_list_find_idx(n,m,num_mol) ][0][a+1]
					H_obj.ExFo[ H_obj.ExFo_starters[m*num_mol+n] + a ] +=  -V_mat( n, m, 0, 0, 0, a+1 )


	# ExFv: ExFv[i,a,b,c]  m<n: <ib|V|ac>;  m>n: <bi|V|ca>
	for m in range(num_mol):
		for n in range(num_mol): 
			if m < n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						for c in range(Nv[n]):
							# H_obj.ExFv[ H_obj.ExFv_starters[m*num_mol+n] + a*Nv[n]**2 + b*Nv[n] + c ] +=  v_mat_list[ v_list_find_idx(m,n,num_mol) ][ b+1 ][ (a+1)*rec[n]+(c+1) ]
							H_obj.ExFv[ H_obj.ExFv_starters[m*num_mol+n] + a*Nv[n]**2 + b*Nv[n] + c ] += V_mat( m, n, 0, b+1, a+1, c+1 )
			if m > n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						for c in range(Nv[n]):
							# H_obj.ExFv[ H_obj.ExFv_starters[m*num_mol+n] + a*Nv[n]**2 + b*Nv[n] + c ] +=  v_mat_list[ v_list_find_idx(n,m,num_mol) ][ b+1 ][ (c+1)*rec[m]+(a+1)]
							H_obj.ExFv[ H_obj.ExFv_starters[m*num_mol+n] + a*Nv[n]**2 + b*Nv[n] + c ] += V_mat( n, m, b+1, 0, c+1, a+1 )


	# ExDx ExDx[i,a,b,j]  <ib|V|aj>   <bi|V|ja>
	for m in range(num_mol):
		for n in range(num_mol):
			if m < n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						# print("ExDx[%d,%d,%d,%d] = %lf" %(0,a,b,0, -v_mat_list[ v_list_find_idx(m,n,num_mol) ][ b+1 ][ a+1 ]))
						# H_obj.ExDx[ H_obj.ExDx_starters[m*num_mol+n] + a*Nv[n] + b ] +=  -v_mat_list[ v_list_find_idx(m,n,num_mol) ][ b+1 ][ a+1 ]
						H_obj.ExDx[ H_obj.ExDx_starters[m*num_mol+n] + a*Nv[n] + b ] +=  -V_mat( m, n, 0, b+1, a+1, 0 )
			elif m > n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						# print("ExDx[%d,%d,%d,%d] = %lf" %(0,a,b,0, -v_mat_list[ v_list_find_idx(n,m,num_mol) ][ b+1 ][ a+1 ]))
						# H_obj.ExDx[ H_obj.ExDx_starters[m*num_mol+n] + a*Nv[n] + b ] +=  -v_mat_list[ v_list_find_idx(n,m,num_mol) ][ b+1 ][ a+1 ]
						H_obj.ExDx[ H_obj.ExDx_starters[m*num_mol+n] + a*Nv[n] + b ] +=  -V_mat( n, m, b+1, 0, 0, a+1 )


	# FoFo: FoFo[i,j,k,l]  <ik|V|jl>   <ki|V|lj>
	for m in range(num_mol):
		for n in range(num_mol):
			if m > n:
				# H_obj.FoFo[ H_obj.FoFo_starters[m*num_mol+n] ] +=  v_mat_list[ v_list_find_idx(n,m,num_mol) ][0][0]
				H_obj.FoFo[ H_obj.FoFo_starters[m*num_mol+n] ] +=  V_mat( n, m, 0, 0, 0, 0 )


	# FvFv: FvFv[a,b,c,d]   <ac|V|bd>   <ca|V|db>
	for m in range(num_mol):
		for n in range(num_mol): 				
			if m > n:
				for a in range(Nv[m]):
					for b in range(Nv[m]):
						for c in range(Nv[n]):
							for d in range(Nv[n]):
								# H_obj.FvFv[ H_obj.FvFvs_starters[m*num_mol+n] +  a*Nv[m]*Nv[n]**2 + b*Nv[n]**2 + c*Nv[n] + d ] += v_mat_list[ v_list_find_idx(n,m,num_mol) ][(c+1)*rec[m] + a+1 ][ (d+1)*rec[m]+ b+1 ]
								H_obj.FvFv[ H_obj.FvFv_starters[m*num_mol+n] +  a*Nv[m]*Nv[n]**2 + b*Nv[n]**2 + c*Nv[n] + d ] += V_mat( n, m, c+1, a+1, d+1, b+1 )


	# FoDx: FoDx[i,j,a, k]   <ia|V|jk>  <ai|V|kj>
	for m in range(num_mol):
		for n in range(num_mol):
			if m < n:
				for a in range(Nv[n]):
					# H_obj.FoDx[ H_obj.FoDx_starters[m*num_mol+n] + a ] += v_mat_list[ v_list_find_idx(m,n,num_mol) ][ a+1 ][ 0 ]
					H_obj.FoDx[ H_obj.FoDx_starters[m*num_mol+n] + a ] += V_mat( m, n, 0, a+1, 0, 0 )
			elif m > n:
				for a in range(Nv[n]):
					# H_obj.FoDx[ H_obj.FoDx_starters[m*num_mol+n] + a ] += v_mat_list[ v_list_find_idx(n,m,num_mol) ][ a+1 ][ 0 ]
					H_obj.FoDx[ H_obj.FoDx_starters[m*num_mol+n] + a ] += V_mat( n, m, a+1, 0, 0, 0 )
				

	# FvDx: FvDx[a,b,c,i]  <ac|V|bi>  <ca|V|ib>
	for m in range(num_mol):
		for n in range(num_mol): 		
			if m < n:		
				for a in range(Nv[m]):
					for b in range(Nv[m]):
						for c in range(Nv[n]):
							# H_obj.FvDx[ H_obj.FvDx_starters[m*num_mol+n] + a*Nv[m]*Nv[n] + b*Nv[n] + c ] += -v_mat_list[ v_list_find_idx(m,n,num_mol) ][ (a+1)*rec[n] + (c+1) ][ b+1 ]
							H_obj.FvDx[ H_obj.FvDx_starters[m*num_mol+n] + a*Nv[m]*Nv[n] + b*Nv[n] + c ] += -V_mat( m, n, a+1, c+1, b+1, 0 )
			elif m > n:
				for a in range(Nv[m]):
					for b in range(Nv[m]):
						for c in range(Nv[n]):
							# H_obj.FvDx[ H_obj.FvDx_starters[m*num_mol+n] + a*Nv[m]*Nv[n] + b*Nv[n] + c ] += -v_mat_list[ v_list_find_idx(n,m,num_mol) ][ (c+1)*rec[m] + (a+1) ][ b+1 ]
							H_obj.FvDx[ H_obj.FvDx_starters[m*num_mol+n] + a*Nv[m]*Nv[n] + b*Nv[n] + c ] += -V_mat( n, m, c+1, a+1, 0, b+1 )


	# FoFv:  FoFv[i,j,a,b]  <ia|V|jb>  <ai|V|bj>
	for m in range(num_mol):
		for n in range(num_mol):
			if m < n:
				for a in range(Nv[n]):
					for b in range(Nv[n]):
						# H_obj.FoFv[ H_obj.FoFv_starters[m*num_mol+n] + a*Nv[n] + b ] +=  -v_mat_list[ v_list_find_idx(m,n,num_mol) ][ a+1 ][ b+1 ]
						H_obj.FoFv[ H_obj.FoFv_starters[m*num_mol+n] + a*Nv[n] + b ] += -V_mat( m, n, 0, a+1, 0, b+1 )
			elif m > n:
				for a in range(Nv[n]):
					for b in range(Nv[n]):
						# H_obj.FoFv[ H_obj.FoFv_starters[m*num_mol+n] + a*Nv[n] + b ] +=  -v_mat_list[ v_list_find_idx(n,m,num_mol) ][ a+1 ][ b+1 ]
						H_obj.FoFv[ H_obj.FoFv_starters[m*num_mol+n] + a*Nv[n] + b ] += -V_mat( n, m, a+1, 0, b+1, 0 )


	# DxDx: DxDx[a,i,b,j]  <ab|V|ij>  <ba|V|ji>
	for m in range(num_mol):
		for n in range(num_mol):
			if m > n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						# H_obj.DxDx[ H_obj.DxDx_starters[m*num_mol+n] + a*Nv[n] + b ] +=  v_mat_list[ v_list_find_idx(n,m,num_mol) ][(b+1)*rec[m] + a+1][0]
						H_obj.DxDx[ H_obj.DxDx_starters[m*num_mol+n] + a*Nv[n] + b ] += V_mat( n, m, b+1, a+1, 0, 0 )


	return H_obj


class cc_fast_obj(object):
	"""operator class holding H, X, T, or Omega components"""
	def __init__(self):
		self.K      = None
		self.Ex     = None
		self.Fo     = None
		self.Fv     = None
		self.Dx     = None
		self.ExEx   = None   
		self.ExFo   = None   
		self.ExFv   = None   
		self.ExDx   = None
		self.FoFo   = None
		self.FvFv   = None
		self.FoFv   = None   
		self.FoDx   = None   
		self.FvDx   = None   
		self.DxDx   = None
		self.ExFoFo = None
		self.ExFvFv = None
		self.ExFoFv = None
		self.ExFoDx = None
		self.ExFvDx = None
		self.ExExDx = None

def print_amps(cc_obj):
	cc_obj_dict = cc_obj.__dict__
	for key in sorted(cc_obj_dict.keys()):
		if '_' not in key and len(key) <5:
			print(key)
			print(cc_obj_dict[key])

def get_fast_op(fast_obj, rec_num_states ):
	No = len(rec_num_states)
	Nv = 0
	for n in rec_num_states:
		Nv += n-1
	fast_obj.K      = np.zeros(1)
	fast_obj.Ex     = np.zeros((No*Nv))
	fast_obj.Fo     = np.zeros((No*No))
	fast_obj.Fv     = np.zeros((Nv*Nv))
	fast_obj.Dx     = np.zeros((Nv*No))
	fast_obj.ExEx   = np.zeros(((No*Nv)**2))
	fast_obj.ExFo   = np.zeros((Nv*No**3))
	fast_obj.ExFv   = np.zeros((No*Nv**3))
	fast_obj.ExDx   = np.zeros(((No*Nv)**2))
	fast_obj.FoFo   = np.zeros((No**4))
	fast_obj.FvFv   = np.zeros((Nv**4))
	fast_obj.FoDx   = np.zeros((Nv*No**3))
	fast_obj.FvDx   = np.zeros((No*Nv**3))
	fast_obj.DxDx   = np.zeros(((No*Nv)**2))
	fast_obj.FoFv   = np.zeros(((No*Nv)**2))
	# fast_obj.ExFoFo = np.zeros((No**5*Nv))
	# fast_obj.ExFvFv = np.zeros((Nv**5*No))
	# fast_obj.ExFoDx = np.zeros(((No*Nv)**2*No**2))
	# fast_obj.ExFvDx = np.zeros(((No*Nv)**2*Nv**2))
	# fast_obj.ExExDx = np.zeros(((No*Nv)**3))

	return fast_obj


def load_values_into_blocks(rec_num_states, all_mol_eigenvalues, v_mat_list):
	H = init_H_op(cc_operator(), rec_num_states)
	H = load_K(H, all_mol_eigenvalues)
	H = load_eigen_values_into_H(H, rec_num_states, all_mol_eigenvalues)
	H = load_v_matrice_into_H(H, rec_num_states, v_mat_list)

	# H_fast = get_fast_op(cc_fast_obj(), rec_num_states)
	# H_fast = load_H( H,  H_fast, rec_num_states, all_mol_eigenvalues, v_mat_list )
	# H_fast = map_main(rec_num_states, H, H_fast)

	# print_amps(H)

	# print("\n\n\n------------------------------------------------------------------------------")
	# print_amps(H_fast)

	# amp_dict = H_fast.__dict__
	# for key in amp_dict.keys():
	# 	if amp_dict[key] != None:
	# 		np.savetxt('local_'+key+'.txt', amp_dict[key])
	# 	else:
	# 		f= open('local_'+key+'.txt', 'w')
	# 		f.write("None")
	# 		f.close()



	# print("\n\n\n------------------------------------------------------------------------------")
	# raise AssertionError
	return H





















'''
if __name__ == "__main__":
	n = 4
	n_states = 4
	rec = [n_states for i in range(n)]
	all_mol_eigenvalues = [ [j+1.0 for j in range(n_states)] for i in range(n)]
	v_mat_list = []
	for i in range(n):
		for j in range(i+1,n):
			v_mat_list += [ np.arange(1, n_states**4+1, 1,np.float64).reshape((n_states**2, n_states**2)) ]

	# print(rec, all_mol_eigenvalues)
	H = load_values_into_blocks(rec, all_mol_eigenvalues, v_mat_list)
	

	for i in range(n):
		print("H.Fo[%d] = %.1f" %(i,H.Fo[i]), end="\t")
		if i % 10 == 0:
			print(' ')
	print(" ")
	for i in range(len(H.Fv)):
		print("H.Fv[%d] = %.1f" %(i,H.Fv[i]), end="\t")
		if i % 10 == 0:
			print(' ')	
	print(" ")
	for i in range(len(H.ExEx)):
		print("H.ExEx[%d] = %.1f" %(i,H.ExEx[i]), end="\t")
		if i % 10 == 0:
			print(' ')	
	print(" ")
	for i in range(len(H.ExFo)):
		print("H.ExFo[%d] = %.1f" %(i,H.ExFo[i]), end="\t")
		if i % 10 == 0:
			print(' ')	
	print(" ")
	for i in range(len(H.ExFv)):
		print("H.ExFv[%d] = %.1f" %(i,H.ExFv[i]), end="\t")
		if i % 10 == 0:
			print(' ')	
	print(" ")
	for i in range(len(H.ExDx)):
		print("H.ExDx[%d] = %.1f" %(i,H.ExDx[i]), end="\t")
		if i % 10 == 0:
			print(' ')	
	print(" ")
	for i in range(len(H.FoFo)):
		print("H.FoFo[%d] = %.1f" %(i,H.FoFo[i]), end="\t")
		if i % 10 == 0:
			print(' ')	
	print(" ")
	for i in range(len(H.FvFv)):
		print("H.FvFv[%d] = %.1f" %(i,H.FvFv[i]), end="\t")
		if i % 10 == 0:
			print(' ')	
	print(" ")
	for i in range(len(H.FoDx)):
		print("H.FoDx[%d] = %.1f" %(i,H.FoDx[i]), end="\t")
		if i % 10 == 0:
			print(' ')	
	print(" ")
	for i in range(len(H.FvDx)):
		print("H.FvDx[%d] = %.1f" %(i,H.FvDx[i]), end="\t")
		if i % 10 == 0:
			print(' ')	
	print(" ")
	for i in range(len(H.FoFv)):
		print("H.FoFv[%d] = %.1f" %(i,H.FoFv[i]), end="\t")
		if i % 10 == 0:
			print(' ')
	print(" ")
'''











