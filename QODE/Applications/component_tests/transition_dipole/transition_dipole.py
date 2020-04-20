#    (C) Copyright 2018 Anthony D. Dutoi
# 
#    This file is part of Qode.
# 
#    Qode is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    Qode is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with Qode.  If not, see <http://www.gnu.org/licenses/>.
#
# Individual Molecular Transition Dipole Matrix Building Module
#
# The Transition Dipole Matrix is computed according to this equation:
# 
#    < M_p" | X_M | M_p' > 
#  = sum_{i=0}^{n_M} L0 \sqrt{ m_0/m_i } \sum_{a=1}^{n_M} U_{ai} \sqrt(1/2) \sqrt{ 1/ \sqrt{ \tilde{\tilde{k}}_a } }
#    x ( \sqrt{n'_a} \delta_{n"_a, n'_a - 1} + \sqrt{ n'_a + 1 } \delta_{n"_a, n'_a + 1} ) 
#    x ( \delta_{1",1}(M) ... \delta_{a"-1,a'-1} \delta_{a"+1,a'+1}(M) ... \delta_{n"_M, n'_M} )
#
#
#
import numpy as np
from math import sqrt
import copy

#
#  SORTING STATES BY ENERGY
#
def sort_by_energy( list_of_index, list_of_energy ):
    # This function sorts the eigen energies and those corresponding
    # oscillator state numbers. 
    # INPUTS :  all states in a list of list , all eigen energies
    # OUTPUTS:  deepcopy of both sorted lists
    sorted_index   =  copy.deepcopy( list_of_index  )
    sorted_energy  =  copy.deepcopy( list_of_energy )
    if len( sorted_index ) != len( sorted_energy ):
        print("SORT_BY_ENERGY FUNCTION GOT INVALID INPUTS!")
        sys.exit(1)
    else:
        total_pairs = len( sorted_index ) - 1 
        for i in range( total_pairs ):
            for j in range( total_pairs - i ):
                #print("compare",j,"=",sorted_energy[j],"and",j+1,"=",sorted_energy[j+1]) 
                if sorted_energy[j] > sorted_energy[j+1]:
                    tmp_E               =  sorted_energy[j] 
                    sorted_energy[j]    =  sorted_energy[j+1]
                    sorted_energy[j+1]  =  tmp_E
                    tmp_i               =  sorted_index[j]
                    sorted_index[j]     =  sorted_index[j+1]
                    sorted_index[j+1]   =  tmp_i   

    return sorted_index, sorted_energy





#
#  Assertation of deltas for pairs of state numbers.
#
def delta_assert( list_of_states1 , list_of_states2 , num_diff_allowed ):
    # One-Dimensional list of states for 
    #    bra:  list_of_states1
    #    ket:  list_of_states2
    # Number of different states Allowed.
    if len(list_of_states1) == len(list_of_states2):
        num_diff = 0
        i = 0                 # index
        while (i < len(list_of_states1) and num_diff <= num_diff_allowed ):
            if list_of_states1[i] != list_of_states2[i]:
                num_diff += 1
            i += 1
        if num_diff == num_diff_allowed:
            return True
        else:
            return False

    else:
        print("Number Of States Not Alignned. Terrible Error!")
        raise TypeError



#
# Get the SINGLE DIFFERENT STATE NUMBER ( must differ by one state or exit )
#
def one_diff_state_index( list_of_states1 , list_of_states2 ):
    num_diff_allowed = 1
    if delta_assert( list_of_states1 , list_of_states2 , num_diff_allowed ):
        return_index = 0
        seek = True
        while return_index < len( list_of_states1 ) and seek:
            if list_of_states1[ return_index ] == list_of_states2[ return_index ]:
                return_index += 1
            else:
                seek = False
        if return_index < len( list_of_states1 ):
            return return_index
        else:
            print("Terrible Error, main function didn't phase out same states ZERO TERMS")
    else:
        print("Terrible Error, main function didn't phase out different state ZERO TERMS")
        sys.exit(1)



#  INTEGRAL ENGINES:
#
# If Harmonic Oscillator Hamiltonian has thus form:
#    H_M = -1/2 d^2 / d \tilde{\tilde{x}}_i^2(M) + 1/2 \tilde{\tilde{k}}_{ii} \tilde{\tilde{x}}_i^2(M)
#
# Its transition dipole integral has this formula:
#    <n'| \tilde{\tilde{x}} |n> = \sqrt{1/2} \sqrt{ 1/ \sqrt{ \tilde{\tilde{k}} } } 
#                                 x ( \sqrt{n} \delta_{n',n-1} + \sqrt{n+1} \delta_{n',n+1} )
#
#
def core_dipole_integral( first_state_num, second_state_num ):
    # integral = <i''|x|i'> =  sqrt( i' ) delta_i'',i'-1 + sqrt( i'+1 ) delta_i'',i'+1 
    if first_state_num == second_state_num - 1:
        return pow( second_state_num , 0.5 ) 
    elif first_state_num == second_state_num + 1:
        return pow( second_state_num + 1  , 0.5 )
    else:
        return 0.0


def dipole_integral_engine( 
                            bra_states      ,
                            ket_states      ,
                            vectors         ,
                            mass_list       ,
                            rotated_k_list  ,
                            unit_length=1.0 , 
                            unit_mass=1.0   , 
                            unit_hbar=1.0      ):
    # states have to be different for ONLY ONE PAIR for each.
    num_diff_allowed = 1 
    if delta_assert( bra_states, ket_states, num_diff_allowed ):
        diff_index  = one_diff_state_index( bra_states  , ket_states  ) 
        num_ho      = len( bra_states )
        integral    = 0.0
        for i in range( num_ho ):
            # sqrt(1/energy_prefactor) counted in this formula...
            integral += unit_length * sqrt( unit_mass / mass_list[i] )              *\
                        vectors[ i, diff_index ]                                    *\
                        sqrt( 1.0 / ( 2.0*sqrt( rotated_k_list[ diff_index ] ) ) )  *\
                        core_dipole_integral( bra_states[ diff_index ], ket_states[ diff_index ] )  
        return integral
    else:
        return 0.0 


#
# Transition Dipole Matrix Builder Main Routine
# 
def transition_dipole_matrix_main(
                                    mol_index        , 
                                    eigen_energy     , 
                                    vectors          ,
                                    mass_list        ,
                                    rotated_k_list   ,
                                    unit_length=1.0  , 
                                    unit_mass=1.0    , 
                                    unit_hbar=1.0          ):
    # vectors come sorted due to the use of numpy.linalg.eigh function.
    # Index and energy are in the order of looping indecies but not in 
    # the order of energy, therefore at first sorting is needed.
    sorted_index1, sorted_energy1 = sort_by_energy( mol_index, eigen_energy ) 
    num_states  = len( mol_index )

    # Filling in off-diagonal elements
    # number of oscillator states that is different must = 1
    dipole_mat = np.zeros((num_states,num_states))

    x_index = 0
    for i in range( num_states ):
        bra_states  = sorted_index1[ i ]
        y_index = 0
        for j in range( num_states ):
            ket_states = sorted_index1[ j ]
            dipole_mat[ x_index, y_index ] += \
                       dipole_integral_engine( 
                                                bra_states     ,  
                                                ket_states     , 
                                                vectors        ,
                                                mass_list      ,
                                                rotated_k_list ,
                                                unit_length    ,  
                                                unit_mass      , 
                                                unit_hbar       )             
            y_index += 1
        x_index += 1 
    return dipole_mat 



def paired_transition_dipole_mat(
                                    first_mol_index        ,
                                    first_energy           ,
                                    first_vectors          ,
                                    first_mass_list        ,
                                    first_rotated_k_list   ,
                                    second_mol_index       ,
                                    second_energy          ,
                                    second_vectors         ,
                                    second_mass_list       ,
                                    second_rotated_k_list  ,
                                    pair_coupling_const    , 
                                    unit_length=1.0        ,  
                                    unit_mass=1.0          , 
                                    unit_hbar=1.0           ):      
    #  Building Paired Transition Dipole Matrix using the following formul a:
    #
    #     < M_p"| <N_q"| V_{MN} |M_p'> |N_q'>   = k_{MN} < M_p"|X_M|N_q'>  <N_q"|X_N|N_q'>
    #  = energy_prefactor * k_{MN} < M_p"|X_M|N_q'>  <N_q"|X_N|N_q'> 
    #
    #
    first_num_states  = len(first_mol_index)
    second_num_states = len(second_mol_index)
    first_dipole_mat  = transition_dipole_matrix_main(
                                                        first_mol_index        ,
                                                        first_energy           ,
                                                        first_vectors          ,
                                                        first_mass_list        ,
                                                        first_rotated_k_list   ,
                                                        unit_length            ,  
                                                        unit_mass              , 
                                                        unit_hbar               )  

    second_dipole_mat  = transition_dipole_matrix_main(
                                                        second_mol_index       ,
                                                        second_energy          ,
                                                        second_vectors         ,
                                                        second_mass_list       ,
                                                        second_rotated_k_list  ,
                                                        unit_length            ,  
                                                        unit_mass              , 
                                                        unit_hbar               )

    paired_dipole_mat = np.zeros((first_num_states*second_num_states, first_num_states*second_num_states))

    x_idx = 0
    for i in range(first_num_states):
        for j in range(second_num_states):
            y_idx = 0
            for k in range(first_num_states):
                for l in range(second_num_states):
                    paired_dipole_mat[x_idx, y_idx] = first_dipole_mat[i,k] * second_dipole_mat[j,l] 
                    y_idx += 1
            x_idx += 1

    # Multiply the pair_coupling_const
    paired_dipole_mat = pair_coupling_const * paired_dipole_mat

    # Multiply Energy Prefactor
    prefactor = unit_hbar ** 2 / ( unit_mass * unit_length **2 )
    if prefactor != 1.0:
        paired_dipole_mat = prefactor * paired_dipole_mat

    return paired_dipole_mat







