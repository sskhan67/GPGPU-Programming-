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
import sys
sys.path.append('..')
import os
from multiprocessing import Pool
import numpy as np
from qode.coupled_oscillators import oscillator_system as osys 
from qode.coupled_oscillators import hamiltonian_class
import driver_recsys

rec_instr = [ 25,25 ]
def get_obj( distance ):
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
    dipole_coupling = -2.0
    k = dipole_coupling / distance**3


    molecule = osys.OscillatorSystem(
               osys.group_of_oscillators(
                           recursive_list=[ mol1, mol2 ]  ,
                           coupling = [[0.0,  k   ],
                                       [0.0,   0.0 ]]        )  )


    recsys_calc = hamiltonian_class.RecSys_HO_hamiltonian(
                            molecule.get_k_mat()              ,
                            molecule.top_level_groups()       ,
                            molecule.top_level_coupling_mat() ,
                            rec_instr                           ) 

    return molecule, recsys_calc


def run( distance ):
    # driver_recsys.driver_main( molecule , recsys_calc )
    molecule_obj, rec_obj = get_obj( distance )
    E_bound = 7.0
    E_a, E_rec = driver_recsys.driver_main( molecule_obj, rec_obj, E_bound )
    name_a   = str( round( distance , 1 ) ) + ".ana.txt"
    name_rec = str( round( distance , 1 ) ) + ".rec.txt"     
    print( name_a )
    print( name_rec ) 

    np.savetxt( name_a,   E_a  )
    np.savetxt( name_rec, E_rec)


#radius = [ 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 15.0,\
#          20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 ]
radius = [15.0]

pool = Pool(25)
pool.map(run,radius)
pool.close()
pool.join()



