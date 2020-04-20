from copy import deepcopy
from .containers_classes import single_amplitudes_class, double_amplitudes_class, amplitudes_class
from .containers_classes import *

def scale(factor, T1):
    if factor == 1.0:
        pass
    else:
        T1_single_cells = T1.single_amplitudes.amplitudes_container
        cell_number = T1.single_amplitudes.cell_number
        for n in range(cell_number):
            T1_single_cells[n] = factor * T1_single_cells[n]

        T1_double_cells = T1.double_amplitudes.amplitudes_container
        for n_2 in range(cell_number):
            for n_3 in range(cell_number):
                for n_1 in range(cell_number):
                    T1_double_cells[n_2, n_3, n_1] = factor * T1_double_cells[n_2, n_3, n_1]


def add_to(T1, T2, scalar):
    T1_single = T1.single_amplitudes
    T1_double = T1.double_amplitudes
    T2_single = T2.single_amplitudes
    T2_double = T2.double_amplitudes
    L = T1.single_amplitudes.cut_off_limit
    for i in range(-L, L+1):
        T1_single.amplitudes_container[i] += scalar * T2_single.amplitudes_container[i]
    for n_1 in range(-L, L+1):
        for n_2 in range(-L, L+1):
            for n_3 in range(-L, L+1):
                T1_double.amplitudes_container[n_2,n_3,n_1] += scalar * T2_double.amplitudes_container[n_2,n_3,n_1]


def dot(T1, T2):
    T1_singles = T1.single_amplitudes
    T2_singles = T2.single_amplitudes
    n_virtuals = T1.single_amplitudes.number_of_virtual_orbitals
    n_occupied = T1.single_amplitudes.number_of_occupied_orbitals
    L = T1.single_amplitudes.cut_off_limit
    m = 0
    result = 0.0
    for n in range(-L, L+1):
        for a in range(n_virtuals):
            for i in range(n_occupied):
                result += T1_singles(a,n,i,m) * T2_singles(a,n,i,m)
    T1_doubles = T1.double_amplitudes
    T2_doubles = T2.double_amplitudes
    for n_1 in range(-L, L+1):
        for n_2 in range(-L, L+1):
            for n_3 in range(-L, L+1):
                for a in range(n_virtuals):
                    for b in range(n_virtuals):
                        for i in range(n_occupied):
                            for j in range(n_occupied):
                                result += T1_doubles(a,n_3,b,n_2,i,m,j,n_1) * T2_doubles(a,n_3,b,n_2,i,m,j,n_1)
    return result


def copy(T2):
    number_of_virtual_orbitals = T2.single_amplitudes.number_of_virtual_orbitals
    number_of_occupied_orbitals = T2.single_amplitudes.number_of_occupied_orbitals
    cut_off_limit = T2.single_amplitudes.cut_off_limit
    amplitudes = deepcopy(T2.single_amplitudes.amplitudes_container)
    t1 = single_amplitudes_class(number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit)
    t1.amplitudes_container = amplitudes
    t2 = double_amplitudes_class(number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit)
    amplitudes2 = deepcopy(T2.double_amplitudes.amplitudes_container)
    t2.amplitudes_container = amplitudes2
    T1 = amplitudes_class(t1, t2)
    return (T1)


def act_on_vec(inv_F, T):
    L = T.single_amplitudes.cut_off_limit
    n_virtuals = T.single_amplitudes.number_of_virtual_orbitals
    n_occupied = T.single_amplitudes.number_of_occupied_orbitals
    T_singles = T.single_amplitudes.amplitudes_container
    energy_diff_singles = inv_F.single_energies.differences_container
    T_doubles = T.double_amplitudes.amplitudes_container
    energy_diff_doubles = inv_F.double_energies.differences_container
    for n in range(-L, L+1):
        for a in range(n_virtuals):
            for i in range(n_occupied):
               T_singles[n, a, i] = T_singles[n, a, i] / energy_diff_singles[n, a, i] 
    for n_1 in range(-L, L+1):
        for n_2 in range(-L, L+1):
            for n_3 in range(-L, L+1):
                for a in range(n_virtuals):
                    for b in range(n_virtuals):
                        for i in range(n_occupied):
                            for j in range(n_occupied):
                                if a==0 and b==1 and i==0 and j==1:
                                    print('*****************************************')
                                    print('a=', a,'b=',b,'i=', i,'j=', j)
                                    print('t=', T_doubles[n_2, n_3, n_1][a, b, i, j], 'e=', energy_diff_doubles[n_2, n_3, n_1][a, b, i, j], 'inside operator')
                                T_doubles[n_2, n_3, n_1][a, b, i, j] = T_doubles[n_2, n_3, n_1][a, b, i, j] / energy_diff_doubles[n_2, n_3, n_1][a, b, i, j]  
                                if a==0 and b==1 and i==0 and j==1:
                                    print('t(updated)=', T_doubles[n_2, n_3, n_1][a, b, i, j])
                                    print('*****************************************')
    print(T_doubles, 'delta T raw')
    return T

