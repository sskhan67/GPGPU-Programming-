from containers_classes import single_amplitudes_class, double_amplitudes_class, amplitudes_class
from containers_classes import one_body_hamiltonian_class, two_body_hamiltonian_class, hamiltonian_class
from containers_classes import omega_class, fill_up_omega, delta_t, fill_up_one_body_hamiltonian
from containers_classes import fill_up_two_body_hamiltonian,fill_up_single_amplitudes, get_omega
from containers_classes import fill_up_double_amplitudes,inv_F_singles, inv_F_doubles, inv_F_class, fill_up_inv_F
from operations_classes import scale, add_to, dot, copy, act_on_vec
import numpy as np


m = 0

# single amplitudes
virtuals = 4
occupied = 2
limit    = 0
'''
t_ia     = single_amplitudes_class(virtuals, occupied, limit)
t_ia     = fill_up_single_amplitudes(t_ia)
t2       = single_amplitudes_class(virtuals, occupied, limit)
t2       = fill_up_single_amplitudes(t2)
tuple2   = (1,3,2)
n = tuple2[0]
i = tuple2[1]
a = tuple2[2]
number3 = n + i**2 + a**3
print(t_ia(a,n,i,m), 'is this the same number as', number3, 'single amplitudes')
'''

# one body hamiltonian
f_pq = one_body_hamiltonian_class(virtuals, occupied, limit)
f_pq = fill_up_one_body_hamiltonian(f_pq)
tuple3 = (0,2,0)
p = tuple3[0]
q = tuple3[1]
n = tuple3[2]
number4 = p + q**3 + n**4
general = virtuals + occupied
M = np.zeros((general,general))
for i in range(general):
    for j in range(general):
        #print('element=', f_pq(i, n, j, m))
        M[i,j] = f_pq(i,n,j,m)
#print(f_pq(p,n,q,m), "is this the same number as", number4, 'one body hamiltonian')

#print(f_pq.matrix_elements_container, "non-zero elements, H core")
#print(M,'is this the abstract matrix (H core) filled up with zeros?')


# two body
V_pqrs = two_body_hamiltonian_class(virtuals, occupied, limit)
V_pqrs = fill_up_two_body_hamiltonian(V_pqrs)
tuple = (0, 0, 0, 0, 0, 0, 0)
n_2 = tuple[4]
n_3 = tuple[5]
n_1 = tuple[6]
p = tuple[0]
q = tuple[1]
r = tuple[2]
s = tuple[3]
element = n_2 + n_3 ** 2 + n_1 ** 3 + p ** 4 + q ** 5 + r ** 6 + s ** 7
g = np.zeros((general, general, general, general))
for i in range(general):
    for j in range(general):
        for k in range(general):
            for l in range(general):
                g[i,j,k,l] = V_pqrs(i, n_2, j, n_3, k, m, l, n_1 )
#print(V_pqrs(p, n_2, q, n_3, r, m, s, n_1), "is this {},{},{},{},  {},{},{}, element = ".format(p, q, r, s, n_2, n_3, n_1),
#      element, 'two body hamiltonian')
#print(V_pqrs.matrix_elements_container[tuple[4], tuple[5], tuple[6], tuple[0], tuple[1], tuple[2], tuple[3]], 'two body matrix element')
print(V_pqrs.matrix_elements_container[0,0,0,1,1], 'slice of tensor')
print(g[1,1], 'slice of filled up tensor' )
print(g[4,4,4,4], 'is this the same number', V_pqrs.matrix_elements_container[0,0,0,2,2,2,2])

"""
# double amplitudes
tuple4 = (2,0,0,0)
a = tuple4[0]
b = tuple4[1]
i = tuple4[2]
j = tuple4[3]
t_abij = double_amplitudes_class(virtuals, occupied, limit)
tt2    = double_amplitudes_class(virtuals, occupied, limit)
t_abij = fill_up_double_amplitudes(t_abij)
tt2    = fill_up_double_amplitudes(tt2)
doubles_element = n_2 + n_3**2 + n_1**3 + a**4 + b**5 + i**6 + j**7
print( t_abij( a, n_2, b, n_3, i, m, j, n_1), "is this the same number as:", doubles_element, 'double amplitudes')

# omega
H  = hamiltonian_class(f_pq, V_pqrs)
T  = amplitudes_class(t_ia, t_abij)
T2 = amplitudes_class(t2, tt2)
print(T, 'object T')
print(T2, 'object T2')
Om = omega_class(T)
Om = fill_up_omega(Om, H, T)
tuple5 = (1,2,7)
a = tuple5[0]
i = tuple5[1]
n = tuple5[2]
number = i + a**2 + n**3
number2 = f_pq(1,1,2,m)
print(Om.singles(a, n, i, m), 'is this number(omega)?', number)

omega = get_omega(H, T)
omega = fill_up_omega(omega, H, T)
print(omega.singles(a, n, i, m), 'is this number?', number)
print(omega.energy, 'is this 5.4')
# operations
print(dot(T,T2), 'result of the dot product')
print(T2.single_amplitudes.amplitudes_container[0], 'is this 1?')
add_to(T,T2, 2)
print(T2.single_amplitudes.amplitudes_container[0], 'is this 1?')
scale(4.5,T2)
print(T2.single_amplitudes.amplitudes_container[0], "is this 4.5 inside T2")
T3 = copy(T)
print(T, "T object 1")
print(T3, 'T object 2')
print(T3.single_amplitudes.amplitudes_container[0], 'is this 4.5 too')

# inv_F
f1 = inv_F_singles(occupied, virtuals, limit)
f2 = inv_F_doubles(occupied, virtuals, limit)
F  = inv_F_class(f1, f2)
F  = fill_up_inv_F(F)
print(F.single_energies.differences_container[0], "is this M_ia = a-i")

print(T.single_amplitudes.amplitudes_container[0], 'contents of T')
T4 = act_on_vec(F, T)
print(T4.single_amplitudes.amplitudes_container[0], '4.5 * differences ')
"""
