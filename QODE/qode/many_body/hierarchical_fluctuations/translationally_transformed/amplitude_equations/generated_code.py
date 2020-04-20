
from tensors import *


L = 1


upper_lim1 = len( t_ai[0][0] )
upper_lim3 = len( t_ai )
upper_lim5 = len( t_ai[0][0] )
upper_lim7 = len( t_ai )
L_0 = [i for i in range(-L, L+1)]
domain_o = L_0


result = 0
for o in domain_o:
    for i in range(upper_lim1):
        L_o = [i for i in range(o-L, o+L+1)]
        domain_p = list( set(L_o) )
        for p in domain_p:
            for a in range(upper_lim3):
                L_p = [i for i in range(p-L, p+L+1)]
                domain_q = list( set(L_p)  & set(L_o) )
                for q in domain_q:
                    for j in range(upper_lim5):
                        L_q = [i for i in range(q-L, q+L+1)]
                        domain_r = list( set(L_p)  & set(L_q)  & set(L_o) )
                        for r in domain_r:
                            for b in range(upper_lim7):
                                result += 0.5*v_pqrs[i][o][j][q][a][p][b][r]*t_ai[a][p][i][o]*t_ai[b][r][j][q]
print(result)

from tensors import *


L = 1


upper_lim1 = len( t_ai )
upper_lim3 = len( t_ai[0][0] )
upper_lim5 = len( t_ai )
upper_lim7 = len( t_ai[0][0] )
L_0 = [i for i in range(-L, L+1)]
domain_o = L_0


result = 0
for o in domain_o:
    for a in range(upper_lim1):
        L_o = [i for i in range(o-L, o+L+1)]
        domain_p = list( set(L_o) )
        for p in domain_p:
            for i in range(upper_lim3):
                L_p = [i for i in range(p-L, p+L+1)]
                domain_q = list( set(L_p)  & set(L_o) )
                for q in domain_q:
                    for b in range(upper_lim5):
                        L_q = [i for i in range(q-L, q+L+1)]
                        domain_r = list( set(L_p)  & set(L_o)  & set(L_q) )
                        for r in domain_r:
                            for j in range(upper_lim7):
                                result += 0.5*v_pqrs[i][p][j][r][a][o][b][q]*t_ai[a][o][i][p]*t_ai[b][q][j][r]
print(result)
