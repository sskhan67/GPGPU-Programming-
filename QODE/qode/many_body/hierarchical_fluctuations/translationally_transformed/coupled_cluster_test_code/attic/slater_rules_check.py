import numpy as np

# load integrals

h = np.load('../integrals_files/core_matrix_He_sto21g.npy')
g = np.load('../integrals_files/eri_tensor_He_sto21g.npy')
s = np.load('../integrals_files/sto21g_1He_overlap_matrix.npy')

# orthogonalizing integrals

N  = (1-s[0,1]**2)**0.5 # normalization factor
x  = -s[0,1]/N
y  = 1/N
c  = np.array([[1, x], [0, y]])
hp = c.T @ h @ c
print(c.T @ s @ c, 'identity?')

gp = np.zeros((g.shape))
for p in range(len(g)):
    for q in range(len(g)):
        for r in range(len(g)):
            for s in range(len(g)):
                for l in range(len(c)):
                    for v in range(len(c)):
                       for w in range(len(c)):
                           for m in range(len(c)):
                               gp[p,q,r,s] += c[l,p]*c[v,q]*c[w,r]*c[m,s]*g[l,v,w,m]
print(gp, 'eri tensor')

 

