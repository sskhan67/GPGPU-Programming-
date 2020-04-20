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
import numpy as np
import numpy.linalg as LA
import qode.coupled_oscillators as osc


def get_analytical_sln( oscillator_system_object, energy_upper_bound, 
                        unit_length, unit_mass, unit_hbar ):
    k_eff, states_configure, E_analytical, vectors = \
                osc.analytical.analytical_main( 
                                         oscillator_system_object.list_of_masses() ,
                                         oscillator_system_object.get_k_mat()      ,
                                         energy_upper_bound                        , 
                                         unit_length,
                                         unit_mass, 
                                         unit_hbar   )

    return np.sort( E_analytical )



def get_two_group_recouple_sln( recouple_hamiltonian_object ):
    dipole_mat        =  recouple_hamiltonian_object.get_pair_dipole_mat()[0]
    energy1, energy2  =  recouple_hamiltonian_object.get_sorted_eigen_energy()

    x = 0
    for i in range(len(energy1)):
        for j in range(len(energy2)):
            dipole_mat[ x ][ x ] += energy1[i] + energy2[j]
            x += 1

    energies, vectors = LA.eigh( dipole_mat, UPLO='U' )
    sorted_energies = np.sort( energies ) 
    return sorted_energies 







def driver_main( oscillator_system_object , recouple_hamiltonian_object, energy_upper_bound ):
    unit_length        = recouple_hamiltonian_object.unit_length
    unit_mass          = recouple_hamiltonian_object.unit_mass
    unit_hbar          = recouple_hamiltonian_object.unit_hbar    



    E_analytical = get_analytical_sln( oscillator_system_object, 
                                       energy_upper_bound, 
                                       unit_length,
                                       unit_mass,
                                       unit_hbar    )
                                                         

    E_recouple   = get_two_group_recouple_sln( recouple_hamiltonian_object )

    return np.sort(E_analytical) , np.sort(E_recouple) 












