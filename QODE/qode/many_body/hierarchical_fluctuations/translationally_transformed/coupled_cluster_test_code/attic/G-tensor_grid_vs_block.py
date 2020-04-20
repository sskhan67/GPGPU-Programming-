import numpy as np
#from energy_from_amplitudes import g_spatial_to_spin_orb, g_spatial_to_spin_orb_phys, H_spatial_to_spin_orb
from integral_transformer   import g_spatial_to_spin_orb_chem_grid, g_spatial_to_spin_orb_phys_grid, H_spatial_to_spin_orb_grid

# this file gives the same result as reference_energy.py

print('###################################################')
electrons = 2
H = np.load('../integrals_files/biorthogonal_semiMO_sto21g_1He_core_matrix.npy')
g = np.load('../integrals_files/biorthogonal_semiMO_sto21g_1He_eri_tensor.npy')
S = np.load('../integrals_files/sto21g_1He_overlap_matrix.npy')

G_phys = g_spatial_to_spin_orb_phys_grid(g)
G      = g_spatial_to_spin_orb_chem_grid(g)


# Density matrix
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


J      = np.tensordot(G,      D, axes=([2,3], [1,0]))
J_phys = np.tensordot(G_phys, D, axes=([3,1], [0,1]))

print(J_phys, 'physiscit notation')
print(J, 'chemist notation')
print('is this J the same?', np.allclose(J, J_phys))
