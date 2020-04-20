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
import numpy as np
import pickle

np.set_printoptions(linewidth=270,threshold=np.nan,precision=2)



def indexprint(np_mat):
	if len(np_mat.shape) != 2: raise TypeError
	v,h = np_mat.shape
	new_np_mat = np.vstack(( np.hstack((np.zeros((1,1)) , np.arange(h).reshape(1,h))) , np.hstack((np.arange(v).reshape(v,1), np_mat)) ))
	print(new_np_mat)

def nonzero_print_mat(np_mat, dim, n_dim):
	if n_dim == 1:
		for x in range(dim):
			if np_mat[x] != 0:
				print('tensor[%d] = %.3e' %(x, np_mat[x]))
	if n_dim == 2:
		for x in range(dim):
			for y in range(dim):
				if np_mat[x,y] != 0:
					print('tensor[%d,%d] = %.3e' %(x, y, np_mat[x,y]))
	if n_dim == 3:
		for x in range(dim):
			for y in range(dim):
				for z in range(dim):
					if np_mat[x,y,z] != 0:
						print('tensor[%d,%d,%d] = %.3e' %(x, y, z, np_mat[x,y,z]))
	if n_dim == 4:
		for x in range(dim):
			for y in range(dim):
				for z in range(dim):
					for w in range(dim):
						if np_mat[x,y,z,w] != 0:
							print('tensor[%d,%d,%d,%d] = %.3e' %(x, y, z, w, np_mat[x,y,z,w]))


def nonzero_print_op(op_list, n, num_states, dim, n_dim):
	for i in range(n):
		for j in range(n):
			if op_list[i][j]:
				for k in range(num_states[i]):
					for l in range(num_states[j]):
						print("op[{}][{}][{}][{}]".format(i,j,k,l) + '-'*40)
						nonzero_print_mat(op_list[i][j][k][l], dim, n_dim)

def printtestcase():
	# print(states[0][2].tolist())
	# print(states[0][10].tolist())
	print(states[1][36].tolist())
	# print(states[1][59].tolist())
	# print(states[2][213].tolist())
	# print(states[2][300].tolist())


if __name__ == '__main__':
	n  = 1
	NN = 1

	states = [ np.load('1e_states.npy'), np.load('2e_states.npy'), np.load('3e_states.npy') ]



	dense_mat = pickle.load( open('dense_mat.p','rb') )

	a  = dense_mat['a']
	c  = dense_mat['c']
	aa = dense_mat['aa']
	cc = dense_mat['cc']
	ca = dense_mat['ca']
	caa = dense_mat['caa']
	cca = dense_mat['cca']
	ccaa = dense_mat['ccaa']

	# print(' ')
	# print('a ' + '='*38 )
	# printtestcase()
	# nonzero_print_op(a,n,[NN,NN,NN], 18, 1)

	# print(' ')
	# print('c ' + '='*38 )
	# printtestcase()
	# nonzero_print_op(c,n,[NN,NN,NN], 18, 1)

	# print(' ')
	# print('aa ' + '='*37 )
	# printtestcase()
	# nonzero_print_op(aa,n,[NN,NN,NN], 18, 2)

	# print(' ')
	# print('cc ' + '='*37 )
	# printtestcase()
	# nonzero_print_op(cc,n,[NN,NN,NN], 18, 2)

	print(' ')
	print('ca ' + '='*37 )
	printtestcase()
	nonzero_print_op(ca,n,[NN,NN,NN], 18, 2)
	print('='*40)
	indexprint(ca[0][0][0][0])

	# print(' ')
	# print('caa ' + '='*36 )
	# printtestcase()
	# nonzero_print_op(caa,n,[NN,NN,NN], 18, 3)

	# print(' ')
	# print('cca ' + '='*36 )
	# printtestcase()
	# nonzero_print_op(cca,n,[NN,NN,NN], 18, 3)

	# print(' ')
	# print('ccaa ' + '='*35 )
	# printtestcase()
	# nonzero_print_op(ccaa,n,[NN,NN,NN], 18, 4)

	# nonzero_print_op(ca,n,[NN,NN,NN], 18, 2)
	# nonzero_print_op(aa,n,[NN,NN,NN], 18, 2)
	# print('='*80)
	# nonzero_print_op(ccaa,n,[NN,NN,NN], 18, 4)








