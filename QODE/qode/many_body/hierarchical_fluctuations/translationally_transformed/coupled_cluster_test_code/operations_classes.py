from copy import deepcopy
from containers_classes import single_amplitudes_class, double_amplitudes_class, amplitudes_class
from containers_classes import *

def scale(factor, T1):
    if factor == 1.0:
        pass
    else:
        T1.single_amplitudes.amplitudes_container = factor * T1.single_amplitudes.amplitudes_container
        T1.double_amplitudes.amplitudes_container = factor * T1.double_amplitudes.amplitudes_container

def add_to(T1, T2, scalar):
    T1_single = T1.single_amplitudes
    T1_double = T1.double_amplitudes
    T2_single = T2.single_amplitudes
    T2_double = T2.double_amplitudes
    T1_single.amplitudes_container += scalar * T2_single.amplitudes_container
    T1_double.amplitudes_container += scalar * T2_double.amplitudes_container


def dot(T1, T2):
    T1_singles = T1.single_amplitudes
    T2_singles = T2.single_amplitudes
    n_virtuals = T1.single_amplitudes.number_of_virtual_orbitals
    n_occupied = T1.single_amplitudes.number_of_occupied_orbitals
    result = 0.0
    for a in range(n_virtuals):
        for i in range(n_occupied):
            result += T1_singles(a,i) * T2_singles(a,i)
    T1_doubles = T1.double_amplitudes
    T2_doubles = T2.double_amplitudes
    for a in range(n_virtuals):
        for b in range(n_virtuals):
            for i in range(n_occupied):
                for j in range(n_occupied):
                    result += T1_doubles(a,b,i,j) * T2_doubles(a,b,i,j)
    return result


def copy(T2):
    number_of_virtual_orbitals  = T2.single_amplitudes.number_of_virtual_orbitals
    number_of_occupied_orbitals = T2.single_amplitudes.number_of_occupied_orbitals
    amplitudes = deepcopy(T2.single_amplitudes.amplitudes_container)
    t1 = single_amplitudes_class(number_of_virtual_orbitals, number_of_occupied_orbitals)
    t1.amplitudes_container = amplitudes
    t2 = double_amplitudes_class(number_of_virtual_orbitals, number_of_occupied_orbitals)
    amplitudes2 = deepcopy(T2.double_amplitudes.amplitudes_container)
    t2.amplitudes_container = amplitudes2
    T1 = amplitudes_class(t1, t2)
    return (T1)


def act_on_vec(inv_F, T):
    n_virtuals = T.single_amplitudes.number_of_virtual_orbitals
    n_occupied = T.single_amplitudes.number_of_occupied_orbitals
    T_singles = T.single_amplitudes.amplitudes_container
    inv_F_singles = inv_F.single_energies.differences_container
    T_doubles = T.double_amplitudes.amplitudes_container
    inv_F_doubles = inv_F.double_energies.differences_container
    for a in range(n_virtuals):
        for i in range(n_occupied):
            T_singles[a, i] = T_singles[a, i] / inv_F_singles[a, i]
    for a in range(n_virtuals):
        for b in range(n_virtuals):  
            for i in range(n_occupied):
                for j in range(n_occupied):
                    if a==0 and b==1 and i==0 and j==1:
                        print('*****************************************')
                        print('a=', a,'b=',b,'i=', i,'j=', j)
                        print('t=', T_doubles[a, b, i, j], 'e=', inv_F_doubles[a, b, i, j], 'inside operator')
                    T_doubles[a, b, i, j] = T_doubles[a, b, i, j]/inv_F_doubles[a, b, i, j]
                    if a==0 and b==1 and i==0 and j==1:
                        print('t(updated)=', T_doubles[a, b, i, j])
                        print('*****************************************')
    print(T_doubles, 'delta T raw')
    return T

