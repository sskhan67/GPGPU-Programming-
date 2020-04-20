import numpy as np
from copy import deepcopy
from .amplitude_equations                 import delta_t1
from .amplitude_equations.T2_equation     import delta_t2
from .amplitude_equations.energy_equation import delta_E
from .integral_transformer                import g_spatial_to_spin_orb_phys_grid, g_spatial_to_spin_orb_chem_grid, H_spatial_to_spin_orb_grid
from .integral_transformer                import get_semiMO_basis, get_G_semiMO_basis



class single_amplitudes_class:
    def __init__(self, number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit):
        cell_number = 2 * cut_off_limit + 1
        self.number_of_virtual_orbitals = number_of_virtual_orbitals
        self.number_of_occupied_orbitals = number_of_occupied_orbitals
        self.cut_off_limit = cut_off_limit
        self.cell_number = cell_number
        amplitudes_container = np.array(
            [np.zeros((number_of_virtual_orbitals, number_of_occupied_orbitals))] * cell_number, dtype=np.float64)
        self.amplitudes_container = amplitudes_container
        self.block_dims = (number_of_virtual_orbitals, number_of_occupied_orbitals, cell_number)
    def __call__(self, a, n, i, m):
        if m!=0: raise AssertionError	# Boris, I'm suspicious of the m!=0 code below (especially the "abs") . . . please check and talk to me about it.
        if m == 0:
            return self.amplitudes_container[n, a, i]
        elif m != 0 and n == 0:
            return self.amplitudes_container[m, a, i]
        elif m != 0 and n != 0:
            return self.amplitudes_container[abs(m-n), a, i]   
    def raw_storage(self):
        return self.amplitudes_container
    def raw_offset(self, n, _zero_):
        if _zero_!=0: raise AssertionError
        d1, d2, B = self.block_dims
        return (n * d1*d2)

class double_amplitudes_class:
    def __init__(self, number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit):
        cell_number = 2*cut_off_limit+1
        self.virtual_orbitals_number  = number_of_virtual_orbitals
        self.occupied_orbitals_number = number_of_occupied_orbitals
        self.cut_off_limit	= cut_off_limit
        self.cell_number   = cell_number
        t_ijab = np.zeros((number_of_virtual_orbitals, number_of_virtual_orbitals, number_of_occupied_orbitals, number_of_occupied_orbitals), dtype=np.float64)
        amplitudes_container = np.array([[[t_ijab] * cell_number] * cell_number] * cell_number)
        self.amplitudes_container = amplitudes_container
        self.block_dims = (number_of_virtual_orbitals, number_of_virtual_orbitals, number_of_occupied_orbitals, number_of_occupied_orbitals, cell_number)
    def __call__(self, a, n_2, b, n_3, i, m, j, n_1):
        if m!=0: raise AssertionError	# Boris, I'm suspicious of the m!=0 code below (especially the "abs") . . . please check and talk to me about it.
        if m == 0:
            return self.amplitudes_container[n_2, n_3, n_1][a, b, i, j]
        elif m != 0 and n_1 == 0:
            n_1p = -n_1
            n_2p =  n_2 - n_1
            n_3p =  n_3 - n_1
            return self.amplitudes_container[n_2p, n_3p, n_1p][a, b, i,j]
        elif m != 0 and n_2 == 0:
            n_1p = n_1 - n_2
            n_2p = -n_2
            n_3p = n_3 - n_2
            return self.amplitudes_container[n_2p, n_3p, n_1p][a, b, i,j]        
        elif m != 0 and n_3 == 0:
            n_1p = n_1 - n_3
            n_2p = n_2 - n_3
            n_3p = -n_3 
            return self.amplitudes_container[n_2p, n_3p, n_1p][a, b, i,j]        
        elif m != 0 and n_1 != 0 and n_2 != 0 and n_3 != 0:
            n_1p = n_1 - m
            n_2p = n_2 - m
            n_3p = n_3 - m
            return self.amplitudes_container[n_2p, n_3p, n_1p][a, b, i, j] 
    def raw_storage(self):
        return self.amplitudes_container
    def raw_offset(self, n_2, n_3, _zero_, n_1):
        if _zero_!=0: raise AssertionError
        d1, d2, d3, d4, B = self.block_dims
        return ((n_2*B*B + n_3*B + n_1) * d1*d2*d3*d4)

class amplitudes_class:
    def __init__(self, single_amplitudes, double_amplitudes):
        self.single_amplitudes = single_amplitudes
        self.double_amplitudes = double_amplitudes

class one_body_hamiltonian_class:
    def __init__(self, number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit):
        cell_number = 2 * cut_off_limit + 1
        general_orbitals = number_of_occupied_orbitals + number_of_virtual_orbitals
        matrix_elements_container = np.array([np.zeros((general_orbitals, general_orbitals))] * cell_number, dtype=np.float64)
        self.virtual_orbitals_number = number_of_virtual_orbitals
        self.occupied_orbitals_number= number_of_occupied_orbitals
        self.matrix_elements_container = matrix_elements_container
        self.cell_number = cell_number
        self.general_orbitals = general_orbitals
        self.cut_off_limit = cut_off_limit
        self.block_dims = (general_orbitals, general_orbitals, cell_number)
    def __call__(self, p, n, q, m):
        if m!=0: raise AssertionError	# Boris, I'm suspicious of the m!=0 code below (especially the "abs") . . . please check and talk to me about it.
        if m == 0:
             return self.matrix_elements_container[n, p, q]
        elif m != 0 and n == 0:
            return self.matrix_elements_container[m, p, q]
        elif m != 0 and n != 0:
            return self.matrix_elements_container[abs(m-n), p, q]  
    def raw_storage(self):
        return self.matrix_elements_container
    def raw_offset(self, n, _zero_):
        if _zero_!=0: raise AssertionError
        d1, d2, B = self.block_dims
        return (n * d1*d2)
             


class two_body_hamiltonian_class:
    def __init__(self, number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit):
        cell_number = 2 * cut_off_limit + 1
        general_orbitals = number_of_occupied_orbitals + number_of_virtual_orbitals
        V_pqrs = np.zeros((general_orbitals, general_orbitals, general_orbitals, general_orbitals), dtype=np.float64)
        matrix_elements_container = np.array([[[V_pqrs] * cell_number] * cell_number] * cell_number)
        self.matrix_elements_container = matrix_elements_container
        self.cell_number = cell_number
        self.occupied_orbitals_number = number_of_occupied_orbitals
        self.virtual_orbitals_number = number_of_virtual_orbitals
        self.general_orbitals = general_orbitals
        self.cut_off_limit = cut_off_limit
        self.block_dims = (general_orbitals, general_orbitals, general_orbitals, general_orbitals, cell_number)
    def __call__(self, p, n_2, q, n_3, r, m, s, n_1):
        if m!=0: raise AssertionError	# Boris, I'm suspicious of the m!=0 code below (especially the "abs") . . . please check and talk to me about it.
        if m == 0:
            n_1p = n_1 
            n_2p = n_2
            n_3p = n_3
            return self.matrix_elements_container[n_1, n_2, n_3][p, q, r, s]
        elif m != 0 and n_1 == 0:
            n_1p = -n_1
            n_2p =  n_2 - n_1
            n_3p =  n_3 - n_1
            return self.matrix_elements_container[n_1p, n_2p, n_3p][p, q, r, s]
        elif m != 0 and n_2 == 0:
            n_1p =  n_1 - n_2
            n_2p = -n_2
            n_3p =  n_3 - n_2
            return self.matrix_elements_container[n_1p, n_2p, n_3p][p, q, r, s]        
        elif m != 0 and n_3 == 0:
            n_1p =  n_1 - n_3
            n_2p =  n_2 - n_3
            n_3p = -n_3
            return self.matrix_elements_container[n_1p, n_2p, n_3p][p, q, r, s]         
        elif m != 0 and n_1 != 0 and n_2 != 0 and n_3 != 0:
            n_1p = n_1 - m
            n_2p = n_2 - m
            n_3p = n_3 - m
            return self.matrix_elements_container[n_1p, n_2p, n_3p][p, q, r, s]        
    def raw_storage(self):
        return self.matrix_elements_container
    def raw_offset(self, n_2, n_3, _zero_, n_1):
        if _zero_!=0: raise AssertionError
        d1, d2, d3, d4, B = self.block_dims
        return ((n_1*B*B + n_2*B + n_3) * d1*d2*d3*d4)



class hamiltonian_class:
    def __init__(self, one_body, two_body):
        self.one_body_hamiltonian = one_body
        self.two_body_hamiltonian = two_body



class omega_class:
    def __init__(self, single_amplitudes, double_amplitudes):
        self.single_amplitudes = single_amplitudes
        self.double_amplitudes = double_amplitudes
        self.energy            = 0.0

def fill_up_omega(omega, H, T):
    L = omega.single_amplitudes.cut_off_limit
    n_occ = omega.single_amplitudes.number_of_occupied_orbitals
    n_vir = omega.single_amplitudes.number_of_virtual_orbitals
    m  = 0  # This corresponds to the position of anonymous cell, UGLY hack but it'll work for now
    f  = H.one_body_hamiltonian
    v  = H.two_body_hamiltonian
    t1 = T.single_amplitudes
    t2 = T.double_amplitudes
    energy = delta_E(L, n_vir, n_occ, f, v, t1, t2)
    omega.energy = energy

    # fill up single amplitudes
    for term in delta_t1.terms:
        term(omega.single_amplitudes, L, n_occ, n_vir, f, v, t1, t2)

    # fill up double amplitudes
    for n_1 in range(-L, L+1):
        for n_2 in range(-L, L+1):
            for n_3 in range(-L, L+1):
                for i in range(n_occ):
                    for j in range(n_occ):
                        for a in range(n_vir):
                            for b in range(n_vir):
                                omega.double_amplitudes.amplitudes_container[n_2, n_3, n_1][a, b, i, j] = delta_t2(i, j, a, b, n_1, n_2, n_3, L, n_vir, n_occ, f, v, t1, t2)
    return omega

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
    elif MOs == 'converged':
        y = 1
        t = 0.77232514 #t =0.072335
        #t = t_transformed(t) # this is for the change of basis 
        c = S[1,0]
        d = S[1,1]
        x = -(b+t*d)/(a+t*c)
        C = np.array([[1, 0, x, 0], [0, 1, 0, x], [t, 0, y, 0], [0, t, 0, y]])
    elif MOs == 'identity':
        C = np.identity(C.shape[0])
    elif MOs == 'HF_reference':
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
    return D

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


def load_integrals(system, orthogonality, base, cut_off_limit, distance, semiMO="none"):  # you need to replace base => basis
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
           if base == '6-31G' or base == '6-31Gmod':
                  if semiMO == "semiMO":
                      h_isolated = [np.load('./integrals_files/{}_1He_core_matrix.npy'.format(base))]
                      g_isolated = [np.load('./integrals_files/{}_1He_eri_tensor.npy'.format(base))]
                      h_dimers = []
                      g_dimers = []
                      g_trimers= []
                      g_trimers_config  =[]
                      g_tetramers       =[]
                      g_tetramers_config=[]
                      h = [h_isolated, h_dimers]
                      g = [g_isolated, g_dimers, g_trimers, g_trimers_config, g_tetramers, g_tetramers_config] # consider using a dictionary for this
                      s = [np.load('./integrals_files/{}_1He_overlap_matrix.npy'.format(base))]
                      for i in range(1, cut_off_limit+1):
                          h_name = './integrals_files/{}_2He_row_0-{}_core_matrix_{}angs.npy'.format(base,i,distance)
                          g_name = './integrals_files/{}_2He_row_0-{}_eri_tensor_{}angs.npy'.format(base,i,distance)
                          h_dimers.append(np.load(h_name))
                          matrix_test = np.load(h_name)
                          g_dimers.append(np.load(g_name))
                          for j in range(i+1,cut_off_limit+1):
                              g_name = './integrals_files/{}_3He_row_0-{}-{}_eri_tensor_{}angs.npy'.format(base,i, j, distance)
                              g_trimers.append(np.load(g_name))
                              g_trimers_config.append((0,i,j))
                              for k in range(j+1,cut_off_limit+1):
                                  g_name = './integrals_files/{}_4He_row_0-{}-{}-{}_eri_tensor_{}angs.npy'.format(base, i, j, k, distance)
                                  g_tetramers.append(np.load(g_name))
                                  g_tetramers_config.append((0,i,j,k))
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
    
def load_parameters(system, orthogonality, basis, cut_off_limit, distance, semiMO="None"):
    if system == 'helium':
        h, g, s, orthogonality = load_integrals(system, orthogonality, basis, cut_off_limit, distance, semiMO)
        if basis == 'sto-21g':
            c = np.array([[0.6138798, 1.2264156], [0.4741148, -1.2869183]])   # obtained from psi4/sto-21g basis
            spatial_orb_energies = [-0.890753,  1.356669]  #HF energies
        if basis == '6-311G':
            c = np.array([[0.2635583, 0.1379505, 1.4438695], [0.4874868, 1.1539942, -1.6987640], [0.3974455, -1.4279369, 0.6962956]]) # obtained form psi4/6311g basis
            spatial_orb_energies = [-0.9168712, 0.8294336, 6.3135464] #HF energies
            #c[:,1]          = c[:,1] + 2*c[:,2]  #artificial non-orthogonality
        if basis == '6-31G':
            if semiMO == 'semiMO':
                #c = np.identity(2)
                c = np.array([[0.5920813, 1.1498181], [0.5135860, -1.1869588]])
                spatial_orb_energies = [-0.9141266, 1.3998593] # HF energies, you miaght need to replace/calculate actual reference energies
            else:
                c = np.array([[0.5920813, 1.1498181],[0.5135860, -1.1869588]])  #obtained from psi4
                spatial_orb_energies = [-0.9141266, 1.3998593] 
        if basis == '6-31Gmod':
            if semiMO == 'semiMO':
                #c = np.identity(2)
                c = np.array([[1.0000000,  0.0000000, 0.0000000, 0.0000000], 
                              [0.0000000,  0.0000000, 0.0000000, 1.0000000], 
                              [0.0000000,  1.0000000, 0.0000000, 0.0000000], 
                              [0.0000000,  0.0000000, 1.0000000, 0.0000000]]) 
                spatial_orb_energies = [-0.3356906, 2.4386836, 2.4386836, 2.4386836] 
    if system == 'neon':
        h, g, s, orthogonality = load_integrals(system, orthogonality, basis, n_atoms, distance, semiMO)
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
        h, g, s, orthogonality = load_integrals(system, orthogonality, basis, n_atoms, distance, semiMO)
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
        h, g, s, orthogonality = load_integrals(system, orthogonality, basis, n_atoms, distance, semiMO)
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

def get_HF_energy_from_fock_matrix(f_pq, parameters): # assumes orthogonal basis
    h, g, s, c, orthogonality, energies = parameters
    pair_elec = int(f_pq.occupied_orbitals_number/2)
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

def get_MP2_correction_spin_orb(v_pqrs, G, spatial_orb_energies ):
    occ = v_pqrs.occupied_orbitals_number
    eps                  = [val for val in spatial_orb_energies for _ in (0, 1)]   # creates the value of the second degenerate spin orbital
    E_MP2 = 0.0
    K = len(eps)
    for a in range(occ):
        for b in range(occ):
            for r in range(occ,K):
                for s in range(occ,K):
                    integrals = (G[a,b,r,s] - G[a,b,s,r])**2
                    energies  =  eps[a] + eps[b] - eps[r] - eps[s]
                    E_MP2    += (1/4.)*integrals/energies
    return(E_MP2)

def fill_up_one_body_hamiltonian(f_pq, parameters):
    h_list, g_list, s, c_HF, orthogonality, energies = parameters
    limit      = f_pq.cut_off_limit
    n_occupied = f_pq.occupied_orbitals_number
    n_orbitals = f_pq.general_orbitals
    f = f_pq.matrix_elements_container
    n_spatial_orbitals = int(f_pq.general_orbitals/2)
    atom_length = int(n_spatial_orbitals/2)
    if limit > 0:
        h_matrices = h_list[0]+h_list[1]
        c_matrices = [np.identity(h_list[0][0].shape[0]), np.identity(h_list[1][0].shape[0])] 
        g_matrices = g_list[0]+g_list[1]
    else:
        h_matrices = h_list[0]
        c_matrices = [np.identity(h_list[0][0].shape[0])]
        g_matrices = g_list[0]
    n_matrices = len(h_matrices)
    for n in range(n_matrices):
        h = h_matrices[n]
        g = g_matrices[n]
        if n <= 1:
           c = c_matrices[n]
        else:
           c = c_matrices[1]
        C = H_spatial_to_spin_orb_grid(c_HF)  #c_HF
        H = H_spatial_to_spin_orb_grid(h)
        G = g_spatial_to_spin_orb_chem_grid(g)
        H_semiMO = get_semiMO_basis(H, C)
        G_semiMO = get_G_semiMO_basis(G, C)
        D_semiMO = get_density_matrix(s[0], n_occupied, 'identity', H, orthogonality)  # check if inputed matrix is used
        J = np.tensordot(G_semiMO, D_semiMO, axes=([2,3], [1,0]))
        K = np.tensordot(G_semiMO, D_semiMO, axes=([1,2], [0,1]))
        Fp = H_semiMO + J - K
        if Fp.shape[0] == n_orbitals:
           pass
        else:
           Fp = Fp[n_orbitals:, :n_orbitals]
        f[ n] = Fp
        f[-n] = Fp.T
    f_pq.matrix_elements_container = f
    return f_pq

def get_configuration_index(m, n1, n2, n3, configurations_key, selection_mode):
    if selection_mode == 'trimer':
        intersection_lenght = 3
    elif selection_mode == 'tetramer':
        intersection_lenght = 4
    unique = np.unique([m, n1, n2, n3])
    unique_scaled = unique + -unique[0]
    for combination in configurations_key:
       intersection = set(unique_scaled) & set(combination)
       if len(intersection) == intersection_lenght and np.array_equal(unique_scaled, combination):
          return configurations_key.index(combination)

def fill_up_two_body_hamiltonian(v_pqrs, parameters):
    h, g, s, c_HF, orthogonality, spatial_orb_energies = parameters
    n_orb_cell = g[0][0].shape[0]
    v = v_pqrs.matrix_elements_container
    n_generals = v_pqrs.general_orbitals
    limit = v_pqrs.cut_off_limit
    m = 0
    for n_1 in range(-limit, limit+1):
        for n_2 in range(-limit, limit+1):
            for n_3 in range(-limit, limit+1):
                unique_indices = list(np.unique([m, n_1, n_2, n_3]))
                n_unique = len(unique_indices)
                if n_unique == 1:
                     g_monomer = g[0][0]
                     g_monomer = get_G_semiMO_basis(g_monomer, c_HF)
                     G_pqrs = g_spatial_to_spin_orb_phys_grid(g_monomer)
                     print("E(MP2/spin)=", get_MP2_correction_spin_orb(v_pqrs, G_pqrs, spatial_orb_energies))
                     G_pqsr = np.moveaxis(G_pqrs, 2,3)
                     v[n_1, n_2, n_3]= G_pqrs - G_pqsr
                if n_unique == 2:
                     n = np.nonzero(unique_indices)[0][0]
                     dimer_index = unique_indices[n]
                     dimer_index = abs(dimer_index) - 1 
                     n_1p = unique_indices.index(n_1)
                     n_2p = unique_indices.index(n_2)
                     n_3p = unique_indices.index(n_3)
                     mp   = unique_indices.index(m)
                     p_beg = n_2p*n_orb_cell
                     p_end = n_2p*n_orb_cell + n_orb_cell
                     q_beg = n_3p*n_orb_cell
                     q_end = n_3p*n_orb_cell + n_orb_cell
                     r_beg =  mp *n_orb_cell
                     r_end =  mp *n_orb_cell + n_orb_cell
                     s_beg = n_1p*n_orb_cell
                     s_end = n_1p*n_orb_cell + n_orb_cell
                     g_dimer = g[1][dimer_index][p_beg:p_end,q_beg:q_end,r_beg:r_end,s_beg:s_end]
                     g_test = g_dimer
                     g_dimer = np.zeros(g_dimer.shape)
                     g_dimer = get_G_semiMO_basis(g_dimer, c_HF)
                     G_pqrs = g_spatial_to_spin_orb_phys_grid(g_dimer)
                     G_pqsr = np.moveaxis(G_pqrs, 2,3)
                     v[n_1, n_2, n_3]= G_pqrs - G_pqsr
                elif n_unique == 3:
                     max_index = max(unique_indices)
                     min_index = min(unique_indices)
                     if abs(max_index - min_index) <= limit:
                         n = get_configuration_index(m, n_1, n_2, n_3, g[3], 'trimer')
                         np_1 = unique_indices.index(n_1) 
                         np_2 = unique_indices.index(n_2) 
                         np_3 = unique_indices.index(n_3) 
                         mp   = unique_indices.index(m) 
                         p_beg = np_2*n_orb_cell
                         p_end = np_2*n_orb_cell + n_orb_cell
                         q_beg = np_3*n_orb_cell
                         q_end = np_3*n_orb_cell + n_orb_cell
                         r_beg =  mp *n_orb_cell
                         r_end =  mp *n_orb_cell + n_orb_cell
                         s_beg = np_1*n_orb_cell
                         s_end = np_1*n_orb_cell + n_orb_cell
                         g_trimer = g[2][n][p_beg:p_end, q_beg:q_end, r_beg:r_end, s_beg:s_end]
                         g_trimer = get_G_semiMO_basis(g_dimer, c_HF)
                         G_pqrs = g_spatial_to_spin_orb_phys_grid(g_trimer)
                         G_pqsr = np.moveaxis(G_pqrs, 2,3)
                         v[n_1, n_2, n_3]= G_pqrs - G_pqsr 
                elif n_unique == 4:
                     max_index = max(unique_indices)
                     min_index = min(unique_indices)
                     if abs(max_index - min_index) <= limit:
                         n = get_configuration_index(m, n_1, n_2, n_3, g[5], 'tetramer')
                         np_1 = unique_indices.index(n_1) 
                         np_2 = unique_indices.index(n_2) 
                         np_3 = unique_indices.index(n_3) 
                         mp   = unique_indices.index(m) 
                         p_beg = np_2*n_orb_cell
                         p_end = np_2*n_orb_cell + n_orb_cell
                         q_beg = np_3*n_orb_cell
                         q_end = np_3*n_orb_cell + n_orb_cell
                         r_beg =  mp *n_orb_cell
                         r_end =  mp *n_orb_cell + n_orb_cell
                         s_beg = np_1*n_orb_cell
                         s_end = np_1*n_orb_cell + n_orb_cell
                         g_tetramer = g[4][n][p_beg:p_end, q_beg:q_end, r_beg:r_end, s_beg:s_end]
                         g_tetramer = get_G_semiMO_basis(g_tetramer, c_HF)
                         G_pqrs = g_spatial_to_spin_orb_phys_grid(g_tetramer)
                         G_pqsr = np.moveaxis(G_pqrs, 2,3)
                         v[n_1, n_2, n_3]= G_pqrs - G_pqsr 
    return v_pqrs


def fill_up_hamiltonian(H, parameters):
    f_pq   = H.one_body_hamiltonian
    v_pqrs = H.two_body_hamiltonian
    f_pq   = fill_up_one_body_hamiltonian(f_pq, parameters)
    v_pqrs = fill_up_two_body_hamiltonian(v_pqrs, parameters)
    return  H



class inv_F_singles:
        # this class contains the energy denominators for single amplitudes
    def __init__(self, number_of_occupied_orbitals, number_of_virtual_orbitals, cut_off_limit):
        self.cut_off_limit = cut_off_limit
        self.occupied_orbitals_number = number_of_occupied_orbitals
        self.virtual_orbitals_number  = number_of_virtual_orbitals
        cell_number = 2*cut_off_limit + 1
        self.cell_number = cell_number
        differences_container = np.array([np.zeros((number_of_virtual_orbitals, number_of_occupied_orbitals))] * cell_number)
        self.differences_container = differences_container
    def __call__(self, a, n, i, m):
        return self.differences_container[n, i, a]

class inv_F_doubles:
        # this class contains the energy denominators for double amplitudes
    def __init__(self, occupied_orbitals_number, virtual_orbitals_number, cut_off_limit):
        self.cut_off_limit = cut_off_limit
        self.occupied_orbitals_number = occupied_orbitals_number
        self.virtual_orbitals_number = virtual_orbitals_number
        cell_number = 2*cut_off_limit + 1
        self.cell_number = cell_number
        e_ijab = np.zeros((virtual_orbitals_number, virtual_orbitals_number, occupied_orbitals_number,
                           occupied_orbitals_number))
        differences_container = np.array([[[e_ijab] * cell_number] * cell_number] * cell_number)
        self.differences_container = differences_container
    def __call__(self, a, n_2, b, n_3, i, m, j, n_1):
        return self.differences_container[n_2, n_3, n_1][p, q, r, s]

class inv_F_class:
    # this class contains the energy denominators classes for singles and doubles
    def __init__(self, singles_energies, doubles_energies):
        self.single_energies = singles_energies
        self.double_energies = doubles_energies
    def __call__(self, a, n, i, m):
        return self.single_energies.differences_container[n, i, a]


def fill_up_inv_F(inv_F, parameters):
    '''
    This class fill up the energy denominators
    '''
    h, g, s, c, orthogonality, spatial_orb_energies = parameters
    #orb_energies = np.load('HF_orb_energies_1He_sto21g.npy')
    orb_energies = [val for val in spatial_orb_energies for _ in (0, 1)]   # creates the value of the second degenerate spin orbital
    singles_container = inv_F.single_energies.differences_container
    L = inv_F.single_energies.cut_off_limit
    n_virtuals = inv_F.single_energies.virtual_orbitals_number
    n_occupied = inv_F.single_energies.occupied_orbitals_number
    occ_ener = orb_energies[:n_occupied]
    vir_ener = orb_energies[n_occupied:]
    for n in range(-L, L+1):
        for a in range(n_virtuals):    # I am using range instead of list because it allows me to calculate dummy energies to test code
            for i in range(n_occupied):
                singles_container[n, a, i] = vir_ener[a] - occ_ener[i]
    doubles_container = inv_F.double_energies.differences_container
    for n_1 in range(-L, L+1):
        for n_2 in range(-L, L+1):
            for n_3 in range(-L, L+1):
                for a in range(n_virtuals):
                    for b in range(n_virtuals):
                        for i in range(n_occupied):
                            for j in range(n_occupied):
                                doubles_container[n_2, n_3, n_1][a, b, i, j] = vir_ener[a]+vir_ener[b]-occ_ener[i]-occ_ener[j]     
    return inv_F



