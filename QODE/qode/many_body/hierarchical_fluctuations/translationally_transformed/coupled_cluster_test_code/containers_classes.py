import numpy as np
from copy                 import deepcopy
from T2_equation          import delta_t2
from T1_equation          import delta_t
from Energy_equation      import delta_E
from integral_transformer import g_spatial_to_spin_orb_phys_grid, g_spatial_to_spin_orb_chem_grid, H_spatial_to_spin_orb_grid

class single_amplitudes_class:
    def __init__(self, number_of_virtual_orbitals, number_of_occupied_orbitals):
        self.number_of_virtual_orbitals  = number_of_virtual_orbitals
        self.number_of_occupied_orbitals = number_of_occupied_orbitals
        amplitudes_container = np.zeros(( number_of_virtual_orbitals, number_of_occupied_orbitals ))
        #t_calculated = 0.77232514     # amplitude obtained from psi4 | t=0.072335 biort_semiMO basis, from Dutoi's coefficient matrix(HF)
        #for i in range(number_of_occupied_orbitals):
        #    amplitudes_container[i,i] = t_calculated 
        self.amplitudes_container = amplitudes_container
    def __call__(self, a, i):
        return self.amplitudes_container[a, i]



class double_amplitudes_class:
    def __init__(self, number_of_virtual_orbitals, number_of_occupied_orbitals):
        self.virtual_orbitals_number  = number_of_virtual_orbitals
        self.occupied_orbitals_number = number_of_occupied_orbitals
        t_ijab = np.zeros( (number_of_virtual_orbitals, number_of_virtual_orbitals, number_of_occupied_orbitals, number_of_occupied_orbitals) )
        amplitudes_container = t_ijab
        self.amplitudes_container = amplitudes_container
    def __call__(self, a, b, i, j):
        #if  a==0 and b==1 and i==0 and j==1:
            #print('t=', self.amplitudes_container[a, b, i, j], 'inside call function')
        return self.amplitudes_container[a, b, i, j]

class amplitudes_class:
    def __init__(self, single_amplitudes, double_amplitudes):
        self.single_amplitudes = single_amplitudes
        self.double_amplitudes = double_amplitudes

class one_body_hamiltonian_class:
    def __init__(self, number_of_virtual_orbitals, number_of_occupied_orbitals):
        general_orbitals = number_of_occupied_orbitals + number_of_virtual_orbitals
        #dimensions = int(general_orbitals / 2)
        matrix_elements_container = np.zeros((general_orbitals, general_orbitals))
        self.number_of_virtual_orbitals  = number_of_virtual_orbitals
        self.number_of_occupied_orbitals = number_of_occupied_orbitals
        self.matrix_elements_container = matrix_elements_container
        self.general_orbitals = general_orbitals
        #self.dimensions = dimensions
    def __call__(self, p, q):
        '''
        dimensions = self.dimensions
        if p % 2 == 0 and q % 2 == 0:
            p = int(p/2)
            q = int(q/2)
        elif p % 2 == 0 and q % 2 == 1:
            return 0.0
        elif p % 2 == 1 and q % 2 == 0:
            return 0.0
        elif p % 2 == 1 and q % 2 == 1:
            p = int( (p-1)/2 )
            q = int( (q-1)/2 )
        '''
        return self.matrix_elements_container[p, q]  # convention I follow upper index(a,p) = columns; lower index(i,q) = rows
             


class two_body_hamiltonian_class:
    def __init__(self, number_of_virtual_orbitals, number_of_occupied_orbitals):  # you need to create a smart retrieving function
        general_orbitals = number_of_occupied_orbitals + number_of_virtual_orbitals
        #dimensions = int(general_orbitals/2)
        matrix_elements_container = np.zeros((general_orbitals, general_orbitals, general_orbitals, general_orbitals))
        #V_pqrs = np.zeros((dimensions, dimensions, dimensions, dimensions))
        #matrix_elements_container = V_pqrs
        self.number_of_virtual_orbitals  = number_of_virtual_orbitals
        self.number_of_occupied_orbitals = number_of_occupied_orbitals
        self.matrix_elements_container = matrix_elements_container
        self.general_orbitals = general_orbitals
        #self.dimensions = dimensions
    def __call__(self, p, q, r, s):
        #dimensions = self.dimensions
        return self.matrix_elements_container[p, q, r, s]
'''
        if p % 2 == 0 and q % 2 == 0 and r % 2 == 0 and s % 2 == 0:   # there might be a more elegant way than to test every case, but it will do for now
            p = int(p/2)
            q = int(q/2)
            r = int(r/2)
            s = int(s/2)
        elif p % 2 == 1 and q % 2 == 0 and r % 2 == 0 and s % 2 == 0:
            return 0.0
        elif p % 2 == 0 and q % 2 == 1 and r % 2 == 0 and s % 2 == 0:
            return 0.0
        elif p % 2 == 0 and q % 2 == 0 and r % 2 == 1 and s % 2 == 0:
            return 0.0
        elif p % 2 == 0 and q % 2 == 0 and r % 2 == 0 and s % 2 == 1:
            return 0.0
        elif p % 2 == 1 and q % 2 == 1 and r % 2 == 0 and s % 2 == 0:
            return 0.0
        elif p % 2 == 1 and q % 2 == 0 and r % 2 == 1 and s % 2 == 0:
            p = int( (p-1)/2 )
            q = int(q/2)
            r = int( (r-1)/2 )
            s = int(s/2)
        elif p % 2 == 1 and q % 2 == 0 and r % 2 == 0 and s % 2 == 1:
            return 0.0
        elif p % 2 == 0 and q % 2 == 1 and r % 2 == 1 and s % 2 == 0:
            return 0.0
        elif p % 2 == 0 and q % 2 == 1 and r % 2 == 0 and s % 2 == 1:
            p = int(p/2)
            q = int( (q-1)/2 )
            r = int( (r-1)/2 )
            s = int(s/2)
        elif p % 2 == 0 and q % 2 == 0 and r % 2 == 1 and s % 2 == 1:
            return 0.0
        elif p % 2 == 1 and q % 2 == 1 and r % 2 == 1 and s % 2 == 0:
            return 0.0
        elif p % 2 == 1 and q % 2 == 1 and r % 2 == 0 and s % 2 == 1:
            return 0.0
        elif p % 2 == 1 and q % 2 == 0 and r % 2 == 1 and s % 2 == 1:
            return 0.0
        elif p % 2 == 0 and q % 2 == 1 and r % 2 == 1 and s % 2 == 1:
            return 0.0
        elif p % 2 == 1 and q % 2 == 1 and r % 2 == 1 and s % 2 == 1:
            p = int((p-1)/2)
            q = int((q-1)/2)
            r = int((r-1)/2)
            s = int((s-1)/2)
        return self.matrix_elements_container[p, q, r, s]
'''



class hamiltonian_class:
    def __init__(self, one_body, two_body):
        self.one_body_hamiltonian              = one_body
        self.two_body_hamiltonian              = two_body  # g_pqrs



class omega_class:
    def __init__(self, single_amplitudes, double_amplitudes):
        self.single_amplitudes = single_amplitudes
        self.double_amplitudes = double_amplitudes
        self.energy            = 0.0

def fill_up_omega(omega, H, T):
    n_occ = omega.single_amplitudes.number_of_occupied_orbitals
    n_vir = omega.single_amplitudes.number_of_virtual_orbitals
    f  = H.one_body_hamiltonian
    v  = H.two_body_hamiltonian
    t1 = T.single_amplitudes
    t2 = T.double_amplitudes
    energy = delta_E(n_vir, n_occ, f, v, t1, t2)
    omega.energy = energy
    for i in range(n_occ):
        for a in range(n_vir):
            omega.single_amplitudes.amplitudes_container[a, i] = delta_t(a, i, n_vir, n_occ, f, v, t1, t2)  # this is the order in generated code
    for i in range(n_occ):
        for j in range(n_occ):
            for a in range(n_vir):
                for b in range(n_vir):
                    if a == 0 and b == 1 and i == 0 and j == 1:
                       print('===========================================')
                       print('a=', a, 'c=', b, 'i=', i, 'k=', j, 'omega numbering')
                    omega.double_amplitudes.amplitudes_container[a, b, i, j] = delta_t2(i, j, a, b, n_vir, n_occ, f, v, t1, t2)
                    if a == 0 and b == 1 and i == 0 and j == 1:
                       print('t(omega)=', omega.double_amplitudes.amplitudes_container[a, b, i, j])
                       print('===========================================')
    return omega



def t_transformed(t):
    S = np.load('./integrals_files/sto21g_1He_overlap_matrix.npy')
    N  = (1-S[0,1]**2)**0.5 # normalization factor
    tp = t*N/(1 + S[0,1]*t)
    return tp


def get_density_matrix(S, electrons, MOs, C, orthogonality):
    a = S[0,0]
    b = S[0,1]
    if MOs == 'reference':
        # Obtaining coefficient matrix
        #        This comes from solving the matrix equation. The MOs have to be physically orthogonal.
        #                 [1, 0].T @ S @ [x, y] = 0
        #        where S is the overlap matrix
        x = 1
        y = -x*a/b
        C = np.array([[1, 0, x, 0], [0, 1, 0, x], [0, 0, y, 0], [0, 0, 0, y]])
    if MOs == 'converged':
        y = 1
        t = 0.77232514 #t =0.072335
        #t = t_transformed(t) # this is for the change of basis 
        c = S[1,0]
        d = S[1,1]
        x = -(b+t*d)/(a+t*c)
        C = np.array([[1, 0, x, 0], [0, 1, 0, x], [t, 0, y, 0], [0, t, 0, y]])
    if MOs == 'identity':
        C = np.identity(4)
    if MOs == 'HF_reference':
        #C = np.array([[0.6138798, 1.2264156], [0.4741148, -1.2869183]])   # HF orbitals/sto-12g
        #C = np.array([[0.2635583, 0.1379505, 1.4438695], [0.4874868, 1.1539942, -1.6987640], [0.3974455, -1.4279369, 0.6962956]]) # obtained form psi4/6311g basis
        C = H_spatial_to_spin_orb_grid(C)
    if orthogonality == 'orthogonal':
        Cocc     = C[:,:electrons]
        D        = Cocc @ Cocc.T
    if orthogonality == 'biorthogonal':
        C_inv    = np.linalg.inv(C)
        Cocc     = C[:,:electrons]
        Cocc_inv = C_inv[:electrons,:]
        D        = Cocc @ Cocc_inv
    #C_inv    = np.linalg.inv(C)
    #Cocc     = C[:,:electrons]
    #Cocc_inv = C_inv[:electrons,:]
    #       D        = Cocc @ Cocc_inv
    #D        = Cocc @ Cocc.T
    return D

def orthogonilize_Hcore(S, h):
    N  = (1-S[0,1]**2)**0.5 # normalization factor
    x  = -S[0,1]/N
    y  = 1/N
    c  = np.array([[1, x], [0, y]])
    hp = c.T @ h @ c
    return(hp, c)


def get_fock_matrix_in_MO(c, f, orthogonality):
    if orthogonality == 'orthogonal':
        fp = c.T @ f @ c
        return fp
    if orthogonality == 'biorthogonal':
        c_inv = np.linalg.inv(c)
        fp = c_inv @ f @ c
        return fp

def get_eri_tensor_in_MO(c, g, orthogonality):
    if orthogonality == 'orthogonal':
        gp = get_eri_tensor_in_ortho_MO(c, g)
        return gp
    if orthogonality == 'biorthogonal':
        gp = get_eri_tensor_in_biortho_MO(c, g)
        return gp

def get_eri_tensor_in_ortho_MO(c, g): 
    g1  = np.tensordot(c, g, axes=([0], [3]))
    g2  = np.tensordot(c, g1,axes=([0], [3]))
    g3  = np.tensordot(c, g2,axes=([0], [3]))
    gp  = np.tensordot(c, g3,axes=([0], [3]))
    return(gp)

def get_eri_tensor_in_biortho_MO(c, g):  
    c_inv = np.linalg.inv(c)
    g1  = np.tensordot(c    , g, axes=([0], [3]))
    g2  = np.tensordot(c_inv, g1,axes=([1], [3]))
    g3  = np.tensordot(c    , g2,axes=([0], [3]))
    gp  = np.tensordot(c_inv, g3,axes=([1], [3]))
    return(gp)


def load_integrals(system, orthogonality, base, semiMO="none"):
    if system == 'helium':
        if orthogonality == 'orthogonal':
           orthogonality_mode = 'orthogonal'
           if base == 'sto-21g':
                  h       = np.load('./integrals_files/core_matrix_He_sto21g.npy')
                  g       = np.load('./integrals_files/eri_tensor_He_sto21g.npy')
                  s       = np.load('./integrals_files/sto21g_1He_overlap_matrix.npy')
           if base == '6-311G':
                  h       = np.load('./integrals_files/6311g_1He_core_matrix.npy')
                  g       = np.load('./integrals_files/6311g_1He_eri_tensor.npy')
                  s       = np.load('./integrals_files/6311g_1He_overlap_matrix.npy')
        if orthogonality == 'biorthogonal':
           orthogonality_mode = 'biorthogonal'
           if base == 'sto-21g':
                  h       = np.load("./integrals_files/biorthogonal_sto21g_1He_core_matrix.npy")
                  g       = np.load("./integrals_files/biorthogonal_sto21g_1He_eri_tensor.npy") 
                  s       = np.load('./integrals_files/sto21g_1He_overlap_matrix.npy')
                  if semiMO == "semiMO":
                      h      = np.load("./integrals_files/biorthogonal_semiMO_sto21g_1He_core_matrix.npy")
                      g      = np.load("./integrals_files/biorthogonal_semiMO_sto21g_1He_eri_tensor.npy")
                      s      = np.load('./integrals_files/sto21g_1He_overlap_matrix.npy')
           if base == '6-311G':
                  h      = np.load("./integrals_files/biorthogonal_6311g_1He_core_matrix.npy")
                  g      = np.load("./integrals_files/biorthogonal_6311g_1He_eri_tensor.npy")
                  s      = np.load('./integrals_files/6311g_1He_overlap_matrix.npy')
    if system == 'neon':
        if orthogonality == 'orthogonal':
            orthogonality_mode = orthogonality
            if base == '6-311G':
                h = np.load('./integrals_files/6311g_1Ne_core_matrix.npy')
                g = np.load('./integrals_files/6311g_1Ne_eri_tensor.npy')
                s = np.load('./integrals_files/6311g_1Ne_overlap_matrix.npy')
            if base == '6-31G':
                h = np.load('./integrals_files/6-31G_1Ne_core_matrix.npy')
                g = np.load('./integrals_files/6-31G_1Ne_eri_tensor.npy')
                s = np.load('./integrals_files/6-31G_1Ne_overlap_matrix.npy')
        if orthogonality == 'biorthogonal':
            orthogonality_mode = orthogonality
            if base == '6-311G':
                h = np.load('./integrals_files/biorthogonal_6311g_1Ne_core_matrix.npy')
                g = np.load('./integrals_files/biorthogonal_6311g_1Ne_eri_tensor.npy')
                s = np.load('./integrals_files/6311g_1Ne_overlap_matrix.npy')
            if base == '6-31G':
                h = np.load('./integrals_files/biorthogonal_6-31G_1Ne_core_matrix.npy')
                g = np.load('./integrals_files/biorthogonal_6-31G_1Ne_eri_tensor.npy')
                s = np.load('./integrals_files/6-31G_1Ne_overlap_matrix.npy')
    if system == 'beryllium':
        if orthogonality == 'orthogonal':
            orthogonality_mode = orthogonality
            if base == 'sto-21g':
                h = np.load('./integrals_files/sto-21g_1Be_core_matrix.npy')
                g = np.load('./integrals_files/sto-21g_1Be_eri_tensor.npy')
                s = np.load('./integrals_files/sto-21g_1Be_overlap_matrix.npy')
            if base == '6-31G':
                h = np.load('./integrals_files/6-31G_1Be_core_matrix.npy')
                g = np.load('./integrals_files/6-31G_1Be_eri_tensor.npy')
                s = np.load('./integrals_files/6-31G_1Be_overlap_matrix.npy')
        if orthogonality == 'biorthogonal':
            orthogonality_mode = orthogonality
            if base == 'sto-21g':
                h = np.load('./integrals_files/biorthogonal_sto-21g_1Be_core_matrix.npy')
                g = np.load('./integrals_files/biorthogonal_sto-21g_1Be_eri_tensor.npy')
                s = np.load('./integrals_files/sto-21g_1Be_overlap_matrix.npy')
            if base == '6-31G':
                h = np.load('./integrals_files/biorthogonal_6-31G_1Be_core_matrix.npy')
                g = np.load('./integrals_files/biorthogonal_6-31G_1Be_eri_tensor.npy')
                s = np.load('./integrals_files/6-31G_1Be_overlap_matrix.npy')
    if system == 'beryllium-helium':
        if orthogonality == 'orthogonal':
            orthogonality_mode = orthogonality
            if base == 'sto-21g':
                h = np.load('./integrals_files/sto-21g_Be-He_core_matrix.npy')
                g = np.load('./integrals_files/sto-21g_Be-He_eri_tensor.npy')
                s = np.load('./integrals_files/sto-21g_Be-He_overlap_matrix.npy')
        if orthogonality == 'biorthogonal':
            orthogonality_mode = orthogonality
            if base == 'sto-21g':
                h = np.load('./integrals_files/biorthogonal_sto-21g_Be-He_core_matrix.npy')
                g = np.load('./integrals_files/biorthogonal_sto-21g_Be-He_eri_tensor.npy')
                s = np.load('./integrals_files/sto-21g_Be-He_overlap_matrix.npy')
    return(h,g,s, orthogonality_mode)
    
def load_parameters(system, orthogonality, basis, semiMO="None"):
    if system == 'helium':
        h, g, s, orthogonality = load_integrals(system, orthogonality, basis, semiMO)
        if basis == 'sto-21g':
            c = np.array([[0.6138798, 1.2264156], [0.4741148, -1.2869183]])   # obtained from psi4/sto-21g basis
            spatial_orb_energies = [-0.890753,  1.356669]  #HF energies
        if basis == '6-311G':
            c = np.array([[0.2635583, 0.1379505, 1.4438695], [0.4874868, 1.1539942, -1.6987640], [0.3974455, -1.4279369, 0.6962956]]) # obtained form psi4/6311g basis
            spatial_orb_energies = [-0.9168712, 0.8294336, 6.3135464] #HF energies
            #c[:,1]          = c[:,1] + 2*c[:,2]  #artificial non-orthogonality
    if system == 'neon':
        h, g, s, orthogonality = load_integrals(system, orthogonality, basis, semiMO)
        if basis == '6-311G':
            c = np.array(
               [[0.5444429, 0.1286608, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1160148, 0.0000000, 0.0000000, 0.0000000, 2.3724806],
                [0.4757234, 0.2173740, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.2331762, 0.0000000, 0.0000000, 0.0000000,-2.4251557],
                [0.0000000, 0.0000000, 0.0000000, 0.3287202, 0.0000000, 0.0000000, 0.0000000, 0.2436075, 0.0000000, 0.0000000, 0.0000000, 1.2448928, 0.0000000],
                [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.3287202, 0.0000000, 0.2436075, 0.0000000, 0.0000000, 0.0000000, 1.2448928, 0.0000000, 0.0000000],
                [0.0000000, 0.0000000, 0.3287202, 0.0000000, 0.0000000, 0.2436075, 0.0000000, 0.0000000, 0.0000000, 1.2448928, 0.0000000, 0.0000000, 0.0000000],
                [0.0062914,-0.6629373, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,-1.4964201, 0.0000000, 0.0000000, 0.0000000, 0.1356036],
                [0.0000000, 0.0000000, 0.0000000, 0.4828994, 0.0000000, 0.0000000, 0.0000000, 0.7581987, 0.0000000, 0.0000000, 0.0000000,-1.3526884, 0.0000000],
                [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.4828994, 0.0000000, 0.7581987, 0.0000000, 0.0000000, 0.0000000,-1.3526884, 0.0000000, 0.0000000],
                [0.0000000, 0.0000000, 0.4828994, 0.0000000, 0.0000000, 0.7581987, 0.0000000, 0.0000000, 0.0000000,-1.3526884, 0.0000000, 0.0000000, 0.0000000],
               [-0.0002833,-0.4645188, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.4545827, 0.0000000, 0.0000000, 0.0000000,-0.0928219],
                [0.0000000, 0.0000000, 0.0000000, 0.4147422, 0.0000000, 0.0000000, 0.0000000,-1.1402132, 0.0000000, 0.0000000, 0.0000000, 0.5131401, 0.0000000],
                [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.4147422, 0.0000000,-1.1402132, 0.0000000, 0.0000000, 0.0000000, 0.5131401, 0.0000000, 0.0000000],
                [0.0000000, 0.0000000, 0.4147422, 0.0000000, 0.0000000,-1.1402132, 0.0000000, 0.0000000, 0.0000000, 0.5131401, 0.0000000, 0.0000000, 0.0000000]])
            spatial_orb_energies =  [-32.7599560,-1.9194359,-0.8415274,-0.8415274,-0.8415274, 1.4096858, 1.4096858, 1.4096858, 1.5905508, 8.0955769, 8.0955769, 8.0955769,86.8730702] 
        if basis == '6-31G':
            c = np.array(
                        [[ 0.9955080, 0.2448237, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1261715],
                         [ 0.0200819,-0.5455531, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,-1.4644450],
                         [ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.6902442, 0.0000000, 0.0000000, 0.9190328, 0.0000000],
                         [ 0.0000000, 0.0000000, 0.6902442, 0.0000000, 0.0000000, 0.0000000, 0.9190328, 0.0000000, 0.0000000],
                         [ 0.0000000, 0.0000000, 0.0000000, 0.6902442, 0.0000000, 0.9190328, 0.0000000, 0.0000000, 0.0000000],
                         [-0.0037833,-0.5485072, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.4352923],
                         [ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.4593170, 0.0000000, 0.0000000,-1.0536063, 0.0000000],
                         [ 0.0000000, 0.0000000, 0.4593170, 0.0000000, 0.0000000, 0.0000000,-1.0536063, 0.0000000, 0.0000000],
                         [ 0.0000000, 0.0000000, 0.0000000, 0.4593170, 0.0000000,-1.0536063, 0.0000000, 0.0000000, 0.0000000]])
            spatial_orb_energies = [-32.7593235, -1.9108192, -0.8307707, -0.8307707, -0.8307707, 1.7557989, 1.7557989, 1.7557989, 1.9702979]
    if system == 'beryllium':
        h, g, s, orthogonality = load_integrals(system, orthogonality, basis, semiMO)
        if basis == 'sto-21g':
            c = np.array([[0.9928982, 0.2938807, 0.0000000, 0.0000000, 0.0000000],
                          [0.0261377,-1.0351471, 0.0000000, 0.0000000, 0.0000000],
                          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0000000],
                          [0.0000000, 0.0000000, 1.0000000, 0.0000000, 0.0000000],
                          [0.0000000, 0.0000000, 0.0000000, 1.0000000, 0.0000000]])
            spatial_orb_energies = [-4.4839921, -0.2540377, 0.2210860, 0.2210860, 0.2210860] 
        if basis == '6-31G':
            c = np.array(
                         [[0.9979958, 0.2226599, 0.0000000, 0.0000000, 0.0000000, 0.0033855, 0.0000000, 0.0000000, 0.0000000],  
                          [0.0161427,-0.2887668, 0.0000000, 0.0000000, 0.0000000, 2.0172448, 0.0000000, 0.0000000, 0.0000000],
                          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.2641521, 0.0000000, 0.0000000, 0.0000000, 1.3194024],
                          [0.0000000, 0.0000000, 0.2641521, 0.0000000, 0.0000000, 0.0000000, 1.3194024, 0.0000000, 0.0000000],
                          [0.0000000, 0.0000000, 0.0000000, 0.2641521, 0.0000000, 0.0000000, 0.0000000, 1.3194024, 0.0000000],
                         [-0.0054832,-0.7609796, 0.0000000, 0.0000000, 0.0000000,-1.8983197, 0.0000000, 0.0000000, 0.0000000],
                          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.8037973, 0.0000000, 0.0000000, 0.0000000,-1.0791241],
                          [0.0000000, 0.0000000, 0.8037973, 0.0000000, 0.0000000, 0.0000000,-1.0791241, 0.0000000, 0.0000000],
                          [0.0000000, 0.0000000, 0.0000000, 0.8037973, 0.0000000, 0.0000000, 0.0000000,-1.0791241, 0.0000000]])
            spatial_orb_energies = [-4.7068905, -0.3012954, 0.0824353,  0.0824353, 0.0824353, 0.4397543, 0.4649310, 0.4649310, 0.4649310]
    if system == 'beryllium-helium':
        h, g, s, orthogonality = load_integrals(system, orthogonality, basis, semiMO)
        if basis == 'sto-21g':
            c = np.array(
                         [[ 0.9928182, 0.0073805, 0.2944963, 0.0000000, 0.0000000, 0.0052679, 0.0261208],
                          [ 0.0264414,-0.0239178,-1.0374317, 0.0000000, 0.0000000,-0.0179976,-0.0980508],
                          [ 0.0006177,-0.0169770, 0.0111304, 0.0000000, 0.0000000,-1.0037911,-0.1588359],
                          [ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0000000, 0.0000000, 0.0000000],
                          [ 0.0000000, 0.0000000, 0.0000000, 1.0000000, 0.0000000, 0.0000000, 0.0000000],
                          [ 0.0011962,-0.6131994, 0.0486580, 0.0000000, 0.0000000, 0.0885526,-1.2304787],
                          [-0.0029316,-0.4706760, 0.0472212, 0.0000000, 0.0000000, 0.0451116, 1.3183654]])
            spatial_orb_energies = [-4.4812064, -0.8976373, -0.2502392, 0.2223353, 0.2223353, 0.2308303, 1.3873506 ]
    return(h, g, s, c, orthogonality, spatial_orb_energies)


def get_HF_energy_from_spatial_integrals_in_MO(f_pq, parameters): # assumes orthogonal basis
    h, g, s, c, orthogonality, energies = parameters
    pair_elec = int(f_pq.number_of_occupied_orbitals/2)
    hp = c.T @ h @ c
    gp = get_eri_tensor_in_ortho_MO(c, g)
    E_0 = 2*hp[0,0] + gp[0,0,0,0]
    return(E_0)

def get_HF_energy_from_spatial_integrals(f_pq, parameters): # assumes orthogonal basis
    h, g, s, c, orthogonality, energies = parameters
    pair_elec = int(f_pq.number_of_occupied_orbitals/2)
    c_occ = c[:,:pair_elec]
    d  = c_occ @ c_occ.T
    j  = np.tensordot(g ,d, axes=([2,3], [1,0]))
    k  = np.tensordot(g ,d, axes=([1,2], [0,1]))
    f = h + 2.0*j - k
    E_0 = np.tensordot(h+f, d, axes=(([1,0],[0,1])))
    return(E_0)

def get_HF_energy_from_fock_matrix(f_pq, parameters): # assumes orthogonal basis
    h, g, s, c, orthogonality, energies = parameters
    pair_elec = int(f_pq.number_of_occupied_orbitals/2)
    c_occ = c[:,:pair_elec]
    # method 1
    d  = c_occ @ c_occ.T
    j  = np.tensordot(g ,d, axes=([2,3], [1,0]))
    k  = np.tensordot(g ,d, axes=([1,2], [0,1]))
    f = h + 2.0*j - k
    # method 2
    hp = c.T @ h @ c
    gp = get_eri_tensor_in_ortho_MO(c, g)
    jp = np.zeros(h.shape)
    for p in range(len(h)):
        for q in range(len(h)):
            for i in range(pair_elec):
                jp[p,q] = gp[p, q, i, i]
    kp = np.zeros(h.shape)
    for p in range(len(h)):
        for q in range(len(h)):
            for i in range(pair_elec):
                kp[p,q] = gp[p, i, i, q]
    fpp = hp + 2.0*jp - kp
    fp = c.T @ f @ c
    #print(np.allclose(fpp, fp,), 'are these the same')
    E_0p = hp[0,0] + fp[0,0]
    return(E_0p)

def get_HF_energy_from_fock_matrix_spin_orb(H, G, F, D, C): # assumes orthogonal basis
    Hp = C.T @ H @ C
    Gp = get_eri_tensor_in_ortho_MO(C, G)
    E_3 = 0.5*np.tensordot(H+F, D, axes=([0,1], [1,0]))
    return E_3

def get_HF_energy_from_fock_matrix_NO_density_matrix(H, F, C):  # assumes orthogonal basis
    Fp = C.T @ F @ C
    Hp = C.T @ H @ C
    E_0 = 0.0
    for i in range(2):
        E_0 += Hp[i,i] + Fp[i,i]
    return E_0/2

def get_HF_energy_from_MO_integrals(f_pq, H, G, C): # assumes orthogonal basis
    Hp = C.T @ H @ C
    Gp = get_eri_tensor_in_ortho_MO(C, G)
    one_body = 0.0
    two_body = 0.0
    for a in range(2): #2 electrons
        one_body += Hp[a,a]
    for a in range(2):
        for b in range(2):
            two_body +=  Gp[a, a, b, b] - Gp[a, b, b, a]
    E_0 = one_body + 0.5*two_body
    return E_0

def get_MP2_correction_spatial_orb(v_pqrs, gp, spatial_orb_energies):
    E_mp2 = 0.0
    pair_elec = int(v_pqrs.number_of_occupied_orbitals/2)
    K = len(gp)
    for a in range(pair_elec):
        for b in range(pair_elec):
            for r in range(pair_elec, K):
                for s in range(pair_elec, K):
                    integrals = gp[a,r,b,s]*( (2.0*gp[r,a,s,b]) - gp[r,b,s,a])
                    energies  =  spatial_orb_energies[a] + spatial_orb_energies[b] - spatial_orb_energies[r] - spatial_orb_energies[s]
                    E_mp2    += integrals/energies
    return(E_mp2)

def get_MP2_correction_spin_orb(v_pqrs, G, spatial_orb_energies ):
    print("--------------------------------")
    #print(G)
    occ = v_pqrs.number_of_occupied_orbitals
    eps                  = [val for val in spatial_orb_energies for _ in (0, 1)]   # creates the value of the second degenerate spin orbital
    E_MP2 = 0.0
    K = len(eps)
    for a in range(occ):
        for b in range(occ):
            for r in range(occ,K):
                for s in range(occ,K):
                    #print(Gp[a,r,b,s], Gp[r,a,s,b], Gp[r,b,s,a])
                    integrals = (G[a,b,r,s] - G[a,b,s,r])**2
                    energies  =  eps[a] + eps[b] - eps[r] - eps[s]
                    E_MP2    += (1/4.)*integrals/energies
    return(E_MP2)

def fill_up_one_body_hamiltonian(f_pq, parameters):
    h, g, s, c, orthogonality, energies = parameters
    #print("E_HF (using g)=", get_HF_energy_from_spatial_integrals_in_MO(f_pq, parameters) )
    #print("E_HF (raw int)=", get_HF_energy_from_spatial_integrals(f_pq, parameters))
    #print("E_HF (using f)=", get_HF_energy_from_fock_matrix(f_pq, parameters))
    C       = H_spatial_to_spin_orb_grid(c)
    H       = H_spatial_to_spin_orb_grid(h)
    G       = g_spatial_to_spin_orb_chem_grid(g)
    n_occupied = f_pq.number_of_occupied_orbitals
    D      = get_density_matrix(s, n_occupied, 'HF_reference', c, orthogonality)
    J      = np.tensordot(G, D, axes=([2,3], [1,0]))
    K      = np.tensordot(G, D, axes=([1,2], [0,1]))
    F = H + J - K
    Fp        = get_fock_matrix_in_MO(C, F, orthogonality)
    #Fp = C.T @ F @ C
    #print("E_HF (using F)=", get_HF_energy_from_fock_matrix_spin_orb(H, G, F, D, C))
    #print("E_HF (old sum)=", get_HF_energy_from_fock_matrix_NO_density_matrix(H, F, C))
    #print("E_HF (MO ints)=", get_HF_energy_from_MO_integrals(f_pq, H, G, C))
    f_pq.matrix_elements_container = Fp

    return f_pq


def fill_up_one_body_hamiltonian_biorthogonal(f_pq, parameters):
    h, g, s, c, orthogonality, energies = parameters
    C         = H_spatial_to_spin_orb_grid(c)
    H         = H_spatial_to_spin_orb_grid(h)
    G         = g_spatial_to_spin_orb_chem_grid(g)
    n_occupied= f_pq.number_of_occupied_orbitals
    D         = get_density_matrix(s, n_occupied, 'HF_reference', c, orthogonality)
    J         = np.tensordot(G, D, axes=([2,3], [1,0]))
    K         = np.tensordot(G, D, axes=([1,2], [0,1]))
    F = H + J - K
    Fp        = get_fock_matrix_in_MO(C, F, orthogonality)
    f_pq.matrix_elements_container = Fp
    return f_pq

def fill_up_two_body_hamiltonian(v_pqrs, parameters):
    h, g, s, c, orthogonality, spatial_orb_energies = parameters
    n_generals      = v_pqrs.general_orbitals
    V               = v_pqrs.matrix_elements_container
    gp              = get_eri_tensor_in_MO(c, g, orthogonality) 
    G               = g_spatial_to_spin_orb_phys_grid( gp )  # result in physicist notation
    Gp              = g_spatial_to_spin_orb_chem_grid( gp )
    print('E(MP2/spat)=', get_MP2_correction_spatial_orb(v_pqrs, gp, spatial_orb_energies))
    print("E(MP2/spin)=", get_MP2_correction_spin_orb(v_pqrs, G, spatial_orb_energies))
    for p in range(n_generals):
        for q in range(n_generals):
            for r in range(n_generals):
                for s in range(n_generals):
                    V[p,  q,  r,  s] = (G[p,q,r,s] - G[p,q,s,r])   # using physicist notation
    return v_pqrs


def fill_up_two_body_hamiltonian_biorthogonal(v_pqrs, parameters):
    h, g, s, c, orthogonality, energies = parameters
    n_generals      = v_pqrs.general_orbitals
    V               = np.zeros(v_pqrs.matrix_elements_container.shape)
    gp              = get_eri_tensor_in_MO(c, g, orthogonality)
    #gp              = get_eri_tensor_in_biortho_MO( c, g) 
    G               = g_spatial_to_spin_orb_phys_grid( gp )  # result in physicist notation
    for p in range(n_generals):
        for q in range(n_generals):
            for r in range(n_generals):
                for s in range(n_generals):
                    V[p,  q,  r,  s] = G[p,q,r,s] - G[p,q,s,r]   # using physicist notation
    v_pqrs.matrix_elements_container = V
    return v_pqrs


def fill_up_hamiltonian(H, parameters):
    f_pq   = H.one_body_hamiltonian
    v_pqrs = H.two_body_hamiltonian
    f_pq   = fill_up_one_body_hamiltonian(f_pq, parameters)
    v_pqrs = fill_up_two_body_hamiltonian(v_pqrs, parameters)
    return  H


def fill_up_single_amplitudes(t_ia, amplitude):
    #amplitude = 0.77232514
    t = t_ia.amplitudes_container
    virtual  = t_ia.number_of_virtual_orbitals
    occupied = t_ia.number_of_occupied_orbitals
    for i in range(occupied):
            t[i, i] = amplitude
    return (t_ia)


def get_omega(H, T):        # this function was supposed to be used inside ComputeOmega, probably should be erased
    omega = omega_class(T)
    omega = fill_up_omega(omega, H, T)
    return omega


def fill_up_double_amplitudes(t_abij):
    t = t_abij.amplitudes_container
    n_occupied = t_abij.occupied_orbitals_number
    n_virtuals = t_abij.virtual_orbitals_number
    for a in range(n_virtuals):
        for b in range(a+1, n_virtuals):
            for i in range(n_occupied):
                for j in range(i+1, n_occupied):
                    element = a**4 + b**5 + i**6 + j**7
                    t[a, b, i, j] = 2
    return t_abij



class inv_F_singles:
    def __init__(self, number_of_occupied_orbitals, number_of_virtual_orbitals):
        self.occupied_orbitals_number = number_of_occupied_orbitals
        self.virtual_orbitals_number  = number_of_virtual_orbitals
        differences_container = np.zeros(( number_of_virtual_orbitals, number_of_occupied_orbitals ))
        self.differences_container = differences_container
    def __call__(self, a, i):
        return self.differences_container[a, i]

class inv_F_doubles:
    def __init__(self, occupied_orbitals_number, virtual_orbitals_number):
        self.occupied_orbitals_number = occupied_orbitals_number
        self.virtual_orbitals_number = virtual_orbitals_number
        e_ijab = np.zeros((virtual_orbitals_number, virtual_orbitals_number, occupied_orbitals_number,
                           occupied_orbitals_number))
        differences_container = e_ijab
        self.differences_container = differences_container
    def __call__(self, a, b, i, j):
        return self.differences_container[p, q, r, s]  # is this line not used? a != p

class inv_F_class:
    def __init__(self, singles_energies, doubles_energies):
        self.single_energies = singles_energies
        self.double_energies = doubles_energies
 
  
def fill_up_inv_F(inv_F, parameters):
    h, g, s, c, orthogonality, spatial_orb_energies = parameters
    #spatial_orb_energies = np.load('./integrals_files/reference_orb_energies_1He_sto21g.npy')
    #spatial_orb_energies = [-0.890753,  1.356669]  #HF energies/sto-21g basis
    #spatial_orb_energies = [-0.9168712, 0.8294336, 6.3135464]
    spin_orb_energies = [val for val in spatial_orb_energies for _ in (0, 1)]   # creates the value of the second degenerate spin orbital
    singles_container = inv_F.single_energies.differences_container
    n_virtuals = inv_F.single_energies.virtual_orbitals_number
    n_occupied = inv_F.single_energies.occupied_orbitals_number
    occ_ener = spin_orb_energies[:n_occupied]
    vir_ener = spin_orb_energies[n_occupied:]
    #vir_ener = [val for val in vir_ener for _ in (0, 1)]   # same^
    #n_occupied = len(occ_ener) # spin orbitals
    #n_virtuals = len(vir_ener) 
    for a in range(n_virtuals):    
        for i in range(n_occupied):
            singles_container[a, i] = vir_ener[a] - occ_ener[i]
    doubles_container = inv_F.double_energies.differences_container
    for a in range(n_virtuals):
        for b in range(n_virtuals): 
            for i in range(n_occupied):
                for j in range(n_occupied): 
                    doubles_container[a, b, i, j] = vir_ener[a]+vir_ener[b]-occ_ener[i]-occ_ener[j]    
    return inv_F



