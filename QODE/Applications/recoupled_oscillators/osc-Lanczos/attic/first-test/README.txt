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
To run the code, use:

	python3.3 input.py

Eigenvalues are saved in 'energy_ana.txt'


Please read input.py as an example.

1) First initialize primitive oscillators as shown in input.py

2) Group oscillators into a "molecule" as shown in input.py
       I set up coupling matrix with a variable r so that 
    I can modify the distances easily


3) Set energy_upper_bound = 11.0 Eh

4) Pass all six arguments to the analytical_solution module



This is the main function in module analytical_solution.py:


def analytical_main(   
                      mass_list,
                      primitive_k_mat,
                      energy_upper_bound,
                      unit_length, 
                      unit_mass, 
                      unit_hbar
                    ):

get the following two arguments by calling:

	mass_list        --->   molecule.list_of_masses()
	primitive_k_mat  --->   molecule.get_k_mat()














