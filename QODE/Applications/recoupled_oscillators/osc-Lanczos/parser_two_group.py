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
from qode.coupled_oscillators import oscillator_system as osys
from qode.coupled_oscillators import hamiltonian_class

def build_recouple_object( data ):
    # data is a dict type.
    print( "INPUT PARAMETERS:" )
    print( "NUM OF OSCILLATORS     =", data.num_oscillator )
    print( "MASS                   =", data.mass )
    print( "FORCE CONSTANT         =", data.force_constant )
    print( "INTERNAL COUPLING      =", data.internal_coupling )
    print( "DISTANCE               =", data.distance )
    print( "RECOUPLE NUM OF STATES =", data.rec_instruct )
    
    oscillators1  = [ osys.OscillatorSystem( osys.primitive_oscillator( m = data.mass[i] , k = data.force_constant[i] ) ) \
                     for i in range(data.num_oscillator) ]

    int_coupling1 = [[ data.internal_coupling for i in range(data.num_oscillator) ] for j in range( data.num_oscillator ) ]
    
    
    
    oscillators2  = [ osys.OscillatorSystem( osys.primitive_oscillator( m = data.mass[i] , k = data.force_constant[i] ) ) \
                     for i in range(data.num_oscillator) ]

    int_coupling2 = [[ data.internal_coupling for i in range(data.num_oscillator) ] for j in range( data.num_oscillator ) ]


    # Remove the Diagonal Part in Internal Coupling Matrices
    for i in range( data.num_oscillator ):
        int_coupling1[i][i] = 0.0
        int_coupling2[i][i] = 0.0


    fragment1 = osys.OscillatorSystem(
                                      osys.group_of_oscillators(
                                                                recursive_list=oscillators1 ,
                                                                coupling = int_coupling1
                                                                )
                                      )

    fragment2 = osys.OscillatorSystem(
                                  osys.group_of_oscillators(
                                                            recursive_list=oscillators2 ,
                                                            coupling = int_coupling2
                                                            )
                                  )
    
    ext_k = -2.0/ (data.distance**3)

    molecule = osys.OscillatorSystem(
                                     osys.group_of_oscillators(
                                                               recursive_list=[ fragment1, fragment2 ]  ,
                                                               coupling = [[0.0, ext_k ],
                                                                           [0.0, 0.0   ]]
                                                               )
                                     )

    rec_num_states = [ data.rec_instruct, data.rec_instruct ]

    rec_obj = hamiltonian_class.RecSys_HO_hamiltonian(
                                                      molecule.get_k_mat()              ,
                                                      molecule.top_level_groups()       ,
                                                      molecule.top_level_coupling_mat() ,
                                                      rec_num_states
                                                      )
    return molecule, rec_obj


