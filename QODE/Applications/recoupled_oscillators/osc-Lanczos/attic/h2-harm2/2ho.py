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
import numpy as np 
from multiprocessing import Pool
import qode.coupled_oscillators.oscillator_system as osys
import qode.coupled_oscillators.analytical        as analytical




unit_length =  1.0
unit_mass   =  1.0
unit_hbar   =  1.0
energy_upper_bound = 10.0



def run(r):
    ho1 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
    ho2 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
    constant = -2.0
    k = constant / r**3
    molecule = osys.OscillatorSystem( 
                   osys.group_of_oscillators( 
                                              recursive_list=[ ho1, ho2 ]  ,
                                              coupling = [[0.0,  k ],
                                                          [0.0,  0.0]]
                                             )
                                 )
    k_eff, states_configure, E_analytical, vectors = \
                analytical.analytical_main( 
                                         molecule.list_of_masses() ,
                                         molecule.get_k_mat()      ,
                                         energy_upper_bound                        , 
                                         unit_length,
                                         unit_mass, 
                                         unit_hbar   )

    return sorted( E_analytical )[0]



def getcurve( r ):
    E_ground = run( r / 0.529177 )
    return  str( r ) +'\t'+  str(E_ground) + '\n' 

E = []

r = 0.5
while r<= 10.0:
    E += getcurve( r )
    r += 0.1

while r<= 100.0:
    E += getcurve( r )
    r += 1.0

energy = ''
for i in E:
    energy += i
print(energy)


