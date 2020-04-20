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



# GROUP ONE
#
ho1 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
ho2 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 0.316 ) )

coupling1 = [[0.0, -0.178],
             [0.0,    0.0]]

mol1 = osys.OscillatorSystem(
        osys.group_of_oscillators(   recursive_list=[ho1,ho2] ,
                                         coupling = coupling1       )  )

#
# GROUP TWO
#
ho3 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
ho4 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 0.316 ) )


coupling2 = [[0.0, -0.178],
             [0.0,    0.0]]


mol2 = osys.OscillatorSystem(
        osys.group_of_oscillators(   recursive_list=[ho3,ho4] ,
                                         coupling = coupling2       )  )

#  INTERMOLECULAR COUPLING
#
distance  = 10.0
dipole_coupling = -2.0
k = dipole_coupling / distance**3


molecule = osys.OscillatorSystem(
           osys.group_of_oscillators(
                       recursive_list=[ mol1, mol2 ]  ,
                       coupling = [[0.0,  k   ],
                                   [0.0,   0.0 ]]        )  )


# SET UP PARAMETERS
energy_upper_bound = 12.0
unit_length  = 1.0
unit_mass    = 1.0
unit_hbar    = 1.0

#  RECOUPLE SYSTEM INSTRUCTIONS
rec_instr = [ 100 , 100 ]
rec_calc = hamiltonian_class.RecSys_HO_hamiltonian(
                                                   molecule.get_k_mat()              ,
                                                   molecule.top_level_groups()       ,
                                                   molecule.top_level_coupling_mat() ,
                                                   rec_instr                           ) 




