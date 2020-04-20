# this script calculates the helium reference energy used in CCS using the sto-21g basis

import numpy as np
from integral_transformer import g_spatial_to_spin_orb_chem_grid, g_spatial_to_spin_orb_phys_grid, H_spatial_to_spin_orb_grid 

# load tensors
##H = np.load('./integrals_files/biorthogonal_semiMO_sto21g_1He_core_matrix.npy')
##g = np.load('./integrals_files/biorthogonal_semiMO_sto21g_1He_eri_tensor.npy')
S = np.load('./integrals_files/sto21g_1He_overlap_matrix.npy')

H = np.load('./integrals_files/biorthogonal_sto21g_1He_core_matrix.npy')
g = np.load('./integrals_files/biorthogonal_sto21g_1He_eri_tensor.npy')


n_spatial  = len(H)
#print(n_spatial, 'n_spatial')
electrons  = 2


## Calculating energy using density matrix (method 1)

# creating density matrix
D = np.zeros((n_spatial, n_spatial))
D[0,0] = 1
D[1,1] = 1
#print(D, "density matrix")

J = np.tensordot(g, D, axes=([2,3], [1,0]))
K = np.tensordot(g, D, axes=([1,2], [0,1]))
F = H + 0.5*(J - K)

ref_energies = np.diag(F)
np.save('reference_orb_energies_1He_sto21g.npy', ref_energies)

E_electric = np.tensordot(F, D, axes=([1,0], [0,1]))
E_total    = E_electric                    # since it is single nucleus; repulsion_energy = 0.0
print('E_total (method 1) =', E_total)


# Calculating energy constructing J and K matrices without density matrix (method 2) 

H_ii = 0.0
V_ij = 0.0
for i in range(electrons):
	H_ii += H[i,i]
	for j in range(electrons):
		V_ij += g[i,i,j,j] - g[i,j,j,i]

energy_2  = H_ii + 0.5*V_ij
E_total_2 = energy_2      # mol.nuclear_repulsion_energy() = 0.0 (single nucleus)
print('E_total (method 2) =', E_total_2 )


'''
# test of for loops
G = g_spatial_to_spin_orb(g)
n_orb = len(G)
n_occ = electrons
V = np.zeros((n_orb, n_orb))
for p in range(n_orb):
    for q in range(n_orb):
        for i in range(n_occ):
            V[p,q] = G[p, q, i, i] - G[p, i, i, q]

print(V, 'V matrix test')
print(J-K, 'V matrix correct')
print('is this true?', np.allclose(V, J-K))
'''


G = g_spatial_to_spin_orb_chem_grid(g)
H = H_spatial_to_spin_orb_grid(H)

# Obtaining density matrix

#        This comes from solving the matrix equation
#                 [1, 0].T @ S @ [x, y] = 0
#        where S is the overlap matrix

x = 1
a = S[0,0]
b = S[0,1]
y = -x*a/b

C = np.array([[1, 0, x, 0], [0, 1, 0, x], [0, 0, y, 0], [0, 0, 0, y]])
#print(C, 'C matrix')
C_inv = np.linalg.inv(C)
Cocc = C[:,:electrons]
Cocc_inv = C_inv[:electrons,:]
D = Cocc @ Cocc_inv


J = np.tensordot(G, D, axes=([2,3], [1,0]))
K = np.tensordot(G, D, axes=([1,2], [0,1]))
F = H + 0.5*(J - K)
E_total3 = np.tensordot(F, D, axes=(([0,1],[1,0])))

print('E total (method 3)', E_total3, 'new density matrix')


# test of physicist notation

G_phys = g_spatial_to_spin_orb_phys_grid(g)
H_phys = H_spatial_to_spin_orb_grid(H)

J_phys = np.tensordot(G_phys, D, axes=([3,1], [0,1]))
K_phys = np.tensordot(G_phys, D, axes=([2,1], [0,1]))

print('is this J the same?', np.allclose(J, J_phys))
print('is this K the same?', np.allclose(K, K_phys))

#print(G_phys, 'physicist notation')
#print(G, 'chemist notation')

print(J, 'chemist notation')
print(J_phys, 'physicist notation')


### HF coefficient matrix test ###

C = np.array([[0.6138798, 1.2264156], [0.4741148, -1.2869183]])
print(C, 'coefficient matrix')
print(C.T @ S @ C, 'identity?')
C = H_spatial_to_spin_orb_grid(C)
S = H_spatial_to_spin_orb_grid(S)
print(C.T @ S @ C, 'identity 2?')

'''
# Gram-Smith method
V_12 = C[:,0] @ S @ C[:,1]
V_11 = C[:,0] @ S @ C[:,0]
m1   = C[:,0]
m2   = C[:,1] - V_12/V_11*C[:,0]
print('---------------------')
print(m1, 'm1')
print(m2, 'm2')
print(V_12/V_11, 'number')
print('---------------------')
m1   = np.reshape(m1, (2,1))
m2   = np.reshape(m2, (2,1))
Cp   = np.concatenate((m1, m2), axis=1)
print(Cp, 'Cp')
'''
