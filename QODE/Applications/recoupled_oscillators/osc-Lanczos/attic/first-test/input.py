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
import qode
import qode.coupled_oscillators.oscillator_system as osys

# The numbers being processed in any mathematical physics code always have some implicit units
# associated with them.  These lines give the values of *those* units in some other base unit
# system.  For example, if the code is using the Bohr radius as a unit of length, but you want
# output in SI units, then unit_length should be set to 5.29e-11 (the Bohr radius in meters).
# Similarly, the values of fundamental constants should be specified here in that system.  It
# looks like the below assumes that hbar is the *working* unit of angular momentum.  Should we 
# lift this assumption later?

# Let's say that the trivial choices below reflect that the working units are the same as the 
# printed units, which are atomic units (the last little bit of specification will be needed below).

unit_length =  1.0
unit_mass   =  1.0
unit_hbar   =  1.0



# The OscillatorSystem class is a container for specifying arbitrarily nested groups of groups
# of osillators, including holding only a single primitive oscillator.  One begins by building
# primitive systems, which is done by constructing the OscillatorSystem class with a special 
# "primitive" argument.

#
#  GROUP ONE
#

# Specify four primitive oscillators

ho1 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 3.0 ) )
ho2 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 4.0 ) )
ho3 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
ho4 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.2 ) )

# These lines define the couplings between four primitive oscillators and then implement
# them by building one group (a "molecule") with those couplings, respecting the order given.
# The form of the Hamiltonian for these coupled oscillators in a group is the following
#
# H = [ SUM(i)  Pi^2 /(2*Mi) + (1/2)Ki*Xi^2 ] + [ SUM(i<j)  Vij*Xi*Xj ]
#   = [ SUM(i)  Pi^2 /(2*Mi) ] + [ (1/2)SUM(i,j) Kij*Xi*Xj ]
#
# where Ki=Kii and Kij=Vij are the numbers entered in the upper triangle of the matrix below.
# The recursive_list argument is so named because the members are themselves expected to be
# oscillator systems, meaning that this class is used to recursively build groups of groups, etc.

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

# repeat two more times ...

#
# GROUP TWO
#

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

#
# GROUP THREE
#

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

#
# TOTAL SYSTEM-WIDE COUPLING
#

# At this point, there are three "molecules" in existence, each of which is a collection
# of four coupled harmonic oscillators.  Now we will arrange these molecules in space.
# This is reflected in the strength of the coupling between the molecules here.  Since there
# is only one spatial parameter below, we must be working in one dimension.  Furthermore, the
# oscillators are each 10 distance units from each other.  Clearly, this arrangement is not
# physically possible!  Nevertheless, it demonstrates the code.
#
# In coupling together two "molecules," each oscillator inside each molecule experiences the same 
# coupling to each oscillator inside the other molecule.  This is reasonable because this represents 
# a dipole coupling which is universal, independent of the force constant or mass of the oscillator
# (assuming each oscillator is mimicing fluctuation of the same charge).
#
# In atomic units, as specified above, dipole coupling is equal to -2/r^3.

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


# Debugging print statements

#'''
print(molecule.list_of_masses())
print("==================")
print(molecule.list_of_force_const())
print("==================")
print(molecule.top_level_coupling_mat())
print("==================")
print(molecule.top_level_groups())
print("==================")
qode.util.printline(molecule.get_k_mat(), "\t")
#'''
