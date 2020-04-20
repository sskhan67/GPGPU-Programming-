from tensors import *

L = 2

a = 0
n = 1
i = 0
m = -1

L_m = [i for i in range(m-L, m+L+1)]
L_n = [i for i in range(n-L, n-L+1)]
upper_lim1 = len(t_ai[0][0])

domain_o = list(set(L_m) & set(L_n))
result = 0
for c in range(upper_lim1):
	for o in  domain_o:
		result += f_pq[a][n][c][o]*t_ai[c][o][i][m]
print(result) 
