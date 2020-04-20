from containers_classes import single_amplitudes_class, double_amplitudes_class, amplitudes_class
from containers_classes import one_body_hamiltonian_class, two_body_hamiltonian_class, hamiltonian_class
from containers_classes import omega_class, fill_up_omega, delta_t, fill_up_one_body_hamiltonian
from containers_classes import fill_up_two_body_hamiltonian,fill_up_single_amplitudes, get_omega
from containers_classes import fill_up_double_amplitudes,inv_F_singles, inv_F_doubles, inv_F_class, fill_up_inv_F
from operations_classes import scale, add_to, dot, copy, act_on_vec
import numpy as np


# initial parameters
electrons = 2
occupied = electrons
loaded_matrix = np.load('core_matrix_He_sto21g.npy')
orbitals = int( 2*len(loaded_matrix) )
virtuals = orbitals - occupied
general = virtuals + occupied

'''
# one body hamiltonian
print('one body hamiltonian test')
f_pq = one_body_hamiltonian_class(virtuals, occupied)
f_pq = fill_up_one_body_hamiltonian(f_pq)
M = np.zeros((general,general))
for i in range(general):
    for j in range(general):
        M[i,j] = f_pq(i,j)
print(M, "abstract matrix")
print("")
print(f_pq.matrix_elements_container, "loaded matrix")
print("------------------------------------")
'''

'''
# two body
print('two body hamiltonian test')
V_pqrs = two_body_hamiltonian_class(virtuals, occupied)
V_pqrs = fill_up_two_body_hamiltonian(V_pqrs)
g = np.zeros((general, general, general, general))
for i in range(general):
    for j in range(general):
        for k in range(general):
            for l in range(general):
                g[i,j,k,l] = V_pqrs(i, j, k, l )
print(V_pqrs.matrix_elements_container[0,0], 'slice of loaded tensor')
print(g[0,0], 'slice of abstract tensor' )
print("--------------------------------------")
'''

'''
# single amplitudes
t_ai     = single_amplitudes_class(virtuals, occupied)
t_ai     = fill_up_single_amplitudes(t_ai)
for i in range(occupied):
    for a in range(virtuals):
        number3 = a + i**3
        print(t_ai(a,i), 'is this the same number as', number3, 'single amplitudes')
print("---------------------------")
'''

'''
# double amplitudes
t_abij = double_amplitudes_class(virtuals, occupied)
t_abij = fill_up_double_amplitudes(t_abij)
for i in range(occupied):
    for j in range(occupied):
        for a in range(virtuals):
            for b in range(virtuals):
                doubles_element = a**4 + b**5 + i**6 + j**7
                print( t_abij( a, b, i, j), "is this the same number as:", doubles_element, 'double amplitudes')
print("----------------------------")
'''

'''
# omega
print('omega test')
g  = np.zeros((general, general, general, general))
M  = np.zeros((general, general))
f_pq   = one_body_hamiltonian_class(virtuals, occupied)
f_pq   = fill_up_one_body_hamiltonian(f_pq)
V_pqrs = two_body_hamiltonian_class(virtuals, occupied)
V_pqrs = fill_up_two_body_hamiltonian(V_pqrs)
H      = hamiltonian_class(f_pq, V_pqrs)
t_ai   = single_amplitudes_class(virtuals, occupied)
t_ai   = fill_up_single_amplitudes(t_ai)
t_abij = double_amplitudes_class(virtuals, occupied)
t_abij = fill_up_double_amplitudes(t_abij)
T      = amplitudes_class(t_ai, t_abij)
Om = omega_class(t_ai, t_abij)
#Om = fill_up_omega(Om, H, T)
for a in range(virtuals):
    for i in range(occupied):
        number = a + i**3
        print(Om.single_amplitudes(a,i), 'is this number(omega)?', number)
print("---------------")
'''


# operations

# set up
t_ai   = single_amplitudes_class(virtuals, occupied)
t_ai   = fill_up_single_amplitudes(t_ai)
t_abij = double_amplitudes_class(virtuals, occupied)
t_abij = fill_up_double_amplitudes(t_abij)
T      = amplitudes_class(t_ai, t_abij)

t1 = single_amplitudes_class(virtuals, occupied)
t2 = double_amplitudes_class(virtuals, occupied)
t1 = fill_up_single_amplitudes(t1)
t2 = fill_up_double_amplitudes(t2)
T2 = amplitudes_class(t1, t2)

invf1 = inv_F_singles(occupied, virtuals)
invf2 = inv_F_doubles(occupied, virtuals)
F  = inv_F_class(invf1, invf2)
F  = fill_up_inv_F(F)


'''
# dot function
print('dot test')
print(dot(T,T2), 'result of the dot product')
print(T.single_amplitudes.amplitudes_container, 'is this 1?')
print("---------------")
'''
'''
# add function
print('add test')
scalar1 = 5
add_to(T,T2, scalar1)
print(T.single_amplitudes.amplitudes_container, 'is this {}?'.format( scalar1+1 ))
print("---------------")
'''
'''
# scale function
print('scale test')
T2 = amplitudes_class(t1, t2)
scalar = 6.6
scale(scalar,T2)
print(T2.single_amplitudes.amplitudes_container, "is this {} inside T2".format(scalar))
print("---------------")
'''

'''
# copy function
print('copy test')
T3 = copy(T)
print(T, "T object 1")
scalar2 = 4.5
scale(scalar2,T3)
print(T3, 'T object 2')
print(T3.single_amplitudes.amplitudes_container, 'is this {} too'.format(scalar2))
print(T.single_amplitudes.amplitudes_container, 'this has to be different from ^')
print("-----------------")
'''


# inv_F function
print('inv_f function')
energies = np.load("HF_orb_energies_1He_sto21g.npy")
invf1 = inv_F_singles(occupied, virtuals)
invf2 = inv_F_doubles(occupied, virtuals)
F  = inv_F_class(invf1, invf2)
F  = fill_up_inv_F(F)
print(F.single_energies.differences_container, "is this M_ai = a-i")
print(F.double_energies.differences_container, "is this M_abij = a+b-i-j")
print(T.single_amplitudes.amplitudes_container, 'contents of T')
print('---------------------------')


# act_on_vec function
print('act on vec')
T4 = act_on_vec(F, T)
print(T4.single_amplitudes.amplitudes_container, 'differences' )


