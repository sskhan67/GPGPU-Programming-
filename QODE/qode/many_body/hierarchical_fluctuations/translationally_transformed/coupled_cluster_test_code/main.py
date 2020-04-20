import numpy as np
from baker_campbell_hausdorff import BakerCampbellHausdorff
from utilities import space_traits
from containers_classes import single_amplitudes_class, double_amplitudes_class, amplitudes_class
from containers_classes import one_body_hamiltonian_class, two_body_hamiltonian_class, hamiltonian_class
from containers_classes import omega_class, fill_up_omega, delta_t, fill_up_one_body_hamiltonian
from containers_classes import fill_up_two_body_hamiltonian,fill_up_single_amplitudes, get_omega, fill_up_hamiltonian
from containers_classes import fill_up_double_amplitudes,inv_F_singles, inv_F_doubles, inv_F_class, fill_up_inv_F
from containers_classes import t_transformed, fill_up_one_body_hamiltonian_biorthogonal, fill_up_two_body_hamiltonian_biorthogonal
from containers_classes import load_parameters
from qode.many_body.coupled_cluster.ccsd import CCSD  
#from qode.util.output   import textlog
from qode.util import parallel, output, textlog, indent
from Energy_equation import delta_E
from T1_equation     import delta_t
from T2_equation     import delta_t2

#import matplotlib.pyplot as plt

# this main is for Helium


electrons = 10
parameters = load_parameters('neon', 'orthogonal', '6-31G')

number_of_occupied_orbitals = electrons
number_of_virtual_orbitals  = int(parameters[0].shape[0]*2 - number_of_occupied_orbitals) 


f_pq   = one_body_hamiltonian_class(number_of_virtual_orbitals, number_of_occupied_orbitals)
v_pqrs = two_body_hamiltonian_class(number_of_virtual_orbitals, number_of_occupied_orbitals)
#v_prqs = two_body_hamiltonian_class(number_of_virtual_orbitals, number_of_occupied_orbitals)

#F   = fill_up_one_body_hamiltonian_biorthogonal(f_pq) 
#V = fill_up_two_body_hamiltonian_biorthogonal(v_pqrs)

H = hamiltonian_class(f_pq, v_pqrs)
H = fill_up_hamiltonian(H, parameters)

'''
print("are these Vs the same?:", np.allclose( v_pqrs.matrix_elements_container, V))
print("are these Fs the same?:", np.allclose( f_pq.matrix_elements_container[3,2], F[3,2]))
print(f_pq.matrix_elements_container, 'orthogonal')
print("-----------------------------")
print(F, 'biorthogonal')
'''


#v_pqrs_prqs = H.two_body_antisymmetrized

t_ai = single_amplitudes_class(number_of_virtual_orbitals, number_of_occupied_orbitals)
t_abij = double_amplitudes_class(number_of_virtual_orbitals, number_of_occupied_orbitals)
T = amplitudes_class(t_ai, t_abij)

f1 = inv_F_singles(number_of_occupied_orbitals, number_of_virtual_orbitals)
f2 = inv_F_doubles(number_of_occupied_orbitals, number_of_virtual_orbitals)
inv_F = inv_F_class(f1, f2)
inv_F = fill_up_inv_F(inv_F, parameters)

out, resources = output(log=textlog(echo=True)), parallel.resources(1)   # this was taken from Yuhong's code

BCH = BakerCampbellHausdorff(resources)
E = CCSD(H, T, inv_F, BCH, space_traits, out.log.sublog(), thresh=1e-6)
print( E )

#energy = delta_E(number_of_virtual_orbitals, number_of_occupied_orbitals, f_pq, v_pqrs, t_ai, t_abij)
#print("energy(calculated amplitudes)=", energy)
#amplitudes = np.zeros((number_of_virtual_orbitals, number_of_occupied_orbitals))

'''
################### graphing eq vs inputs (amplitudes) ########################

#def t_transformed(t):
#    S = np.load('./integrals_files/sto21g_1He_overlap_matrix.npy')
#    N  = (1-S[0,1]**2)**0.5 # normailizatin factor
#    tp = t*N/(1 + S[0,1]*t)
#    return tp

limits = (-3, 3) 
#domain_x = [i/10 for i in range( int(limits[0]*10), int(limits[1]*10))  ]
#domain_x.append(0.77232514)
#domain_x.sort()
value = 0.77232514
value = t_transformed(value)
print('value=', value)
domain_x = [value]  # single value input
total_sum = []
term1 = []
term2 = []
term3 = []
term4 = []
term5 = []
term6 = []
term7 = []
term8 = []
term9 = []
term10= []
term11= []
term12= []
term13= []
term14= []
#print('domain x =', domain_x)
for value in domain_x:
#    for i in range(number_of_occupied_orbitals):
        t_ai = fill_up_single_amplitudes(t_ai, value)
        resulting_tuple = delta_t(0, 0, number_of_virtual_orbitals, number_of_occupied_orbitals, f_pq, v_pqrs, t_ai, t_abij)
        print('resulting of T1 =', resulting_tuple)
        #t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, result = resulting_tuple
        #total_sum.append(result)
        #term1.append(t1)
        #term2.append(t2)
        #term3.append(t3)
        #term4.append(t4)
        #term5.append(t5)
        #term6.append(t6)
        #term7.append(t7)
        #term8.append(t8)
        #term9.append(t9)
        #term10.append(t10)
        #term11.append(t11)
        #term12.append(t12)
        #term13.append(t13)
        #term14.append(t14)
#print(total_sum, 'amplitudes (using calculated amplitude)')

data = [domain_x, term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11, term12, term13, term14, total_sum]
#np.save('/integrals_files/data.npy', data)
'''
