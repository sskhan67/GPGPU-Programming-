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
from qode.coupled_oscillators import oscillator_system as osys
from qode.coupled_oscillators import hamiltonian_class
from qode.coupled_oscillators import two_body_hamil_to_three_body
from qode.coupled_oscillators import nbody_hamiltonian_func


ho1 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 3.0 ) )
ho2 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 4.0 ) )
ho3 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
ho4 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.2 ) )


coupling1 = [[0.0, -0.178, -0.178, -0.178],
             [0.0,    0.0, -0.178, -0.178],
             [0.0,    0.0,    0.0, -0.178],
             [0.0,    0.0,    0.0,   0.0]]

mol1 = osys.OscillatorSystem(
            osys.group_of_oscillators(
                                     recursive_list=[ho1,ho2,ho3,ho4],
                                     coupling = coupling1
                                      )
                             )




ho5 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 3.2 ) )
ho6 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 4.4 ) )
ho7 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.3 ) )
ho8 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.1 ) )

coupling2 = [[0.0, -0.178, -0.178, -0.178],
             [0.0,    0.0, -0.178, -0.178],
             [0.0,    0.0,    0.0, -0.178],
             [0.0,    0.0,    0.0,   0.0]]

mol2= osys.OscillatorSystem(
            osys.group_of_oscillators(
                                     recursive_list=[ho5,ho6,ho7,ho8],
                                     coupling = coupling2
                                      )
                             )




ho9  = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 4.4 ) )
ho10 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 5.3 ) )
ho11 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.4 ) )
ho12 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.2 ) )

coupling3 = [[0.0, -0.178, -0.178, -0.178],
             [0.0,    0.0, -0.178, -0.178],
             [0.0,    0.0,    0.0, -0.178],
             [0.0,    0.0,    0.0,   0.0]]

mol3 = osys.OscillatorSystem(
            osys.group_of_oscillators(
                                     recursive_list=[ho9,ho10,ho11,ho12],
                                     coupling = coupling3
                                      )
                             )



constant = -2.0
r = 10.0
k = constant / r**3

molecule = osys.OscillatorSystem(
                   osys.group_of_oscillators(
                                              recursive_list=[ mol1, mol2, mol3 ]  ,
                                              coupling = [[0.0,  k,   k  ],
                                                          [0.0,  0.0, k  ],
                                                          [0.0,  0.0, 0.0]]
                                             )
                                 )

rec_instr = [ 4, 4, 4]

rec_calc = hamiltonian_class.RecSys_HO_hamiltonian( \
                            molecule.get_k_mat()              ,
                            molecule.top_level_groups()       ,
                            molecule.top_level_coupling_mat() ,
                            rec_instr                           ) 



explicit_H = two_body_hamil_to_three_body.two_to_three_main( rec_calc )
#print(explicit_H)


eigen_list = rec_calc.get_sorted_eigen_energy()
coupling_mat_list = rec_calc.get_pair_dipole_mat()

lda = len(explicit_H)

implicit_H = [ [ 0.0 for i in range(lda) ] for j in range(lda) ] 

err_sum = 0.0

x = 0
while x < lda:
    y = 0 
    while y < lda:
        implicit_H[x][y] = nbody_hamiltonian_func.get_hamiltonian_value( eigen_list , coupling_mat_list , x , y )
        error = explicit_H[x][y] - implicit_H[x][y]
        if abs(error) > 1e-10:
            print("[",x,y,"] =",error)
        err_sum += error**2
        y += 1
    x += 1 

print( "sum of error =", err_sum )







