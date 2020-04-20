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
###########################################################
#                      HO parameters                      #
###########################################################

#  Initial inputs for the hamiltonian

states_per_oscillator = 2
mass                  = [1,1,1]
k_OHconstant          = [1,1,1]
k_coupling            = [0.1,0.1,0.1]
hbar                  = 1.0
arbitrary_mass        = 1


# Calculate useful parameters
#
#n_oscillators = len(mass)
#n_states = states_per_oscillator ** n_oscillators
#CI_matrix = np.zeros((n_states,n_states))


