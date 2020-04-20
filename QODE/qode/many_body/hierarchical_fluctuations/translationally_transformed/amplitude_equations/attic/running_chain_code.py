
from tensors import *


L = 1
result = 0

a = 0
n = 0
i = 0
m = 0


upper_lim1 = len( t_ai )
L_n = [i for i in range(n-L, n+L+1)]
L_m = [i for i in range(m-L, m+L+1)]
domain_o = list( set(L_m)  & set(L_n) )


for o in domain_o:
    for c in range(upper_lim1):
        result += 1*f_pq[a][n][c][o]*t_ai[c][o][i][m]





upper_lim1 = len( t_abij[0][0] )
upper_lim3 = len( t_abij )
upper_lim5 = len( t_abij[0][0][0][0] )
upper_lim7 = len( t_abij[0][0][0][0][0][0] )
L_0 = [i for i in range(-L, L+1)]
domain_o = L_0


for o in domain_o:
    for b in range(upper_lim1):
        L_o = [i for i in range(o-L, o+L+1)]
        domain_p = list( set(L_o) )
        for p in domain_p:
            for a in range(upper_lim3):
                L_p = [i for i in range(p-L, p+L+1)]
                domain_q = list( set(L_o)  & set(L_p) )
                for q in domain_q:
                    for i in range(upper_lim5):
                        L_q = [i for i in range(q-L, q+L+1)]
                        domain_r = list( set(L_o)  & set(L_q)  & set(L_p) )
                        for r in domain_r:
                            for j in range(upper_lim7):
                                result += 0.25*v_pqrs[i][q][j][r][a][p][b][o]*t_abij[a][p][b][o][i][q][j][r]
print(result)
