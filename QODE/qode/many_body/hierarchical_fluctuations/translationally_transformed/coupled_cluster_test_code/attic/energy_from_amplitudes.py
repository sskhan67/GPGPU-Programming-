# This files calculates the HF energy from biorthogonal matrices and a coefficient matrix constructed from CC amplitudes
import numpy as np
from containers_classes   import fill_up_one_body_hamiltonian, get_density_matrix
from containers_classes   import one_body_hamiltonian_class
from integral_transformer import H_spatial_to_spin_orb_grid, g_spatial_to_spin_orb_chem_grid, one_electron_to_spin_block, two_electron_to_spin_phys_block, two_electron_to_spin_chem_block
from integral_transformer import g_spatial_to_spin_orb_phys_grid


# load tensors

H = np.load('./integrals_files/biorthogonal_semiMO_sto21g_1He_core_matrix.npy')  # SEMI-MO
g = np.load('./integrals_files/biorthogonal_semiMO_sto21g_1He_eri_tensor.npy')   # SEMO-MO
#H = np.load('./integrals_files/biorthogonal_sto21g_1He_core_matrix.npy')
#g = np.load('./integrals_files/biorthogonal_sto21g_1He_eri_tensor.npy')


electrons = 2

# Calculating energy using cofficient matrix created from amplitudes
# using chemist and grid notation
 
C1     = np.array([[1.0, 0.0, 0.0314949, 0.0],[0.0, 1.0, 0.0, 0.0314949],[0.072335,  0.0, 1.0, 0.0],[0.0, 0.072335,  0.0, 1.0]])    # Re-ordered Dutoi's matrix/grid ordering 
C1_inv = np.linalg.inv(C1)
D1     = C1[:,:electrons] @ C1_inv[:electrons,:]
#print(C, 'coefficient matrix, chemist/grid ordering ')

H1 = H_spatial_to_spin_orb_grid(H)
g1 = g_spatial_to_spin_orb_chem_grid(g)
J1 = np.tensordot(g1, D1, axes=([2,3],[1,0]))   # chemist notation
K1 = np.tensordot(g1, D1, axes=([1,2],[0,1]))
#F1 = H1 + 0.5*(J1-K1)
#E1 = np.tensordot(F1, D1, axes=(([0,1],[1,0])))
F1 = H1 + J1 - K1
E1 = (1/2.) * np.tensordot( H1+F1, D1, axes=(([1,0],[0,1])) )

#print("Helium HF/sto-21g energy =", E1)
print("Energy using chem/grid       notation =", E1)
print(F1, 'Fock matrix grid ordering')
#print(Hp, 'core matrix grid ordering')
#print(Jp, 'culoumb matris grid ordering')
#print(Kp, 'exchange matrix grid ordering')

# using physicist notation and block notation

C2 = np.array([ [1, 0.072335, 0, 0], [0, 0, 1, 0.072335], [0.0314949, 1, 0, 0], [0, 0, 0.0314949, 1] ]).T # Dutoi's block ordering
C2_inv = np.linalg.inv(C2)
D2 = C2[:,:electrons] @ C2_inv[:electrons,:]


H2 = one_electron_to_spin_block(H)
g2 = two_electron_to_spin_phys_block(g)

J2 = np.tensordot(g2, D2, axes=([3,1],[0,1]))   # physicist notation
K2 = np.tensordot(g2, D2, axes=([2,1],[0,1]))
F2 = H2 + J2 - K2
E2 = (1/2.) * np.tensordot( H2+F2, D2, axes=(([1,0],[0,1])) )
print('Energy using physicist/block notation =', E2, "are these energies equal?:", np.allclose(E2, E1))
#print(Fpp, 'block ordering')


print('-------------------------------------------------')

# using chemist notation and block notation

C3 = np.array([ [1, 0.072335, 0, 0], [0, 0, 1, 0.072335], [0.0314949, 1, 0, 0], [0, 0, 0.0314949, 1] ]).T # Dutoi's block ordering
C3_inv = np.linalg.inv(C3)
D3 = C3[:,:electrons] @ C3_inv[:electrons,:]

H3 = one_electron_to_spin_block(H)
g3 = two_electron_to_spin_chem_block(g)

J3 = np.tensordot(g3, D3, axes=([2,3],[1,0])) # chemist notation
K3 = np.tensordot(g3, D3, axes=([1,2],[0,1]))
F3 = H3 + J3 - K3
E3 = (1/2.) * np.tensordot( H3+F3, D3, axes=(([1,0],[0,1])) )
print('Energy using chemist/block   notation =', E3, "are these energies equal?:", np.allclose(E3, E1))
#print(F3, 'block ordering')

# using physicist notation and grid notation
 
C4     = np.array([[1.0, 0.0, 0.0314949, 0.0],[0.0, 1.0, 0.0, 0.0314949],[0.072335,  0.0, 1.0, 0.0],[0.0, 0.072335,  0.0, 1.0]])    # Re-ordered Dutoi's matrix/grid ordering 
C4_inv = np.linalg.inv(C4)
D4 = C4[:,:electrons] @ C4_inv[:electrons,:]
#print(C, 'coefficient matrix physicist/grid ordering')

H4 = H_spatial_to_spin_orb_grid(H)
g4 = g_spatial_to_spin_orb_phys_grid(g)

J4 = np.tensordot(g4, D4, axes=([3,1],[0,1])) # physicist notation
K4 = np.tensordot(g4, D4, axes=([2,1],[0,1]))
F4 = H4 + J4 - K4
E4 = (1/2.) * np.tensordot( H4+F4, D4, axes=(([1,0],[0,1])) )
print('Energy using physicist/grid  notation =', E4, "are these energies equal?:", np.allclose(E4, E1))

#print(F4, 'fock matrix, grid notation')
#print(C4, 'coefficient matrix')
#print(H4, 'core matrix grid ordering')
#print(D4, 'Density matrix/grid ordering')
#print(g4, 'V tensor physicisti notation')
#print(J4, 'culoumb matrix grid ordering')
#print(K4, 'exchange matrix grid ordering')


# test class container
print('---------------------------------')

n_occupied = electrons
n_virtuals = len(H1) - n_occupied
fock = one_body_hamiltonian_class(n_virtuals, n_occupied)
fock = fill_up_one_body_hamiltonian(fock)
n_general = fock.general_orbitals
Fock_mat = np.zeros((n_general, n_general))
for p in range(n_general):
    for q in range(n_general):
        Fock_mat[p, q] = fock(p,q)
print(Fock_mat, 'Fock operator')

#E_total2 = np.tensordot(Fock_mat, Dp, axes=(([0,1],[1,0])))

#print("Helium HF/sto-21g energy =", E_total2, 'test of class')


