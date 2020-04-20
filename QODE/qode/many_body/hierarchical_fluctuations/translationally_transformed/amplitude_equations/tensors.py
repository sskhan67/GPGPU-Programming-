import numpy as np

# This file only use is for testing CC equations. 


#numerical values of tensors
#f_pq   = [[[[1]*d]*d]*d]*d
#v_pqrs = [[[[[[[[ 1 ] *d]*d]*d]*d]*d]*d]*d]*d
#t_ai   = [[[[1]*d]*d]*d]*d
#t_abij = [[[[[[[[ 1 ] *d]*d]*d]*d]*d]*d]*d]*d

n_occupied = 2
n_virtuals = 2
n_generals = n_occupied + n_virtuals
limit      = 2

def f_pq(p, n, q, m):
	#f = [[[[1]*d]*d]*d]*d
	f = np.ones((n_generals, n_generals, n_generals, n_generals))
	#print('a=', a, 'm=', m, 'i=', i, 'n=', n)
	return f[p][n][q][m]


def t_ai(a, m, i, n):
	#t = [[[[1]*d]*d]*d]*d
	t = np.ones((n_generals, n_generals, n_generals, n_generals))
	return t[a][m][i][n]


def v_pqrs(p,n_2,q,n_3,r,m,s,n_1):
	#v = [[[[[[[[ 1 ] *d]*d]*d]*d]*d]*d]*d]*d
	v = np.ones((n_generals, n_generals, n_generals, n_generals, n_generals, n_generals, n_generals, n_generals))
	return v[p][n_2][q][n_3][r][m][s][n_1]


def t_abij(a,n_2,b,n_3,i,m,j,n_1):
	#t = [[[[[[[[ 1 ] *d]*d]*d]*d]*d]*d]*d]*d
	t = np.ones((n_generals, n_generals, n_generals, n_generals, n_generals, n_generals, n_generals, n_generals))
	return t[a][n_2][b][n_3][i][m][j][n_1]
