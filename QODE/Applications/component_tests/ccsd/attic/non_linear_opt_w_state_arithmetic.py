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
from copy import deepcopy
import numpy as         np
from qode.fermion_field import occ_strings
from qode.fermion_field.state import state, dot, resolvent
from hamiltonian import hamiltonian, SD_excitations

# Energy calculated from Q-Chem 4.3:
# SCF energy in the final basis set = -2.8551604262
# CCSD total energy                 = -2.87178917
#
# PSI4 CISD energy                  = -2.8717891505016



occ_orbs = list(reversed([0,4]))
vrt_orbs = list(reversed([1,2,3,5,6,7]))

cisd_strings = occ_strings.CISD(occ_orbs,vrt_orbs)
num_config = len(cisd_strings)
reference = state(cisd_strings,0)
T = SD_excitations(reference, reference, cisd_strings, occ_orbs, vrt_orbs)

h_mat = np.load('data/h_mat.npy')
V_mat = np.load('data/V_mat.npy')
H = hamiltonian(h_mat, V_mat, cisd_strings, occ_orbs, vrt_orbs)

D = deepcopy(H)
Eo = dot(reference,H(reference))
D.C0 -= (Eo - 1e-10)
R = resolvent(D)

debug = True

while True:
	accumulate = state(cisd_strings)
	accumulate.increment(H(reference),-1)
	intermediate = state(cisd_strings)	
	intermediate.increment(H(T(reference)))
	intermediate.increment(T(H(reference)),-1)
	#
	E = Eo + dot(reference,intermediate)
	print("E =", E)
	if debug:
		psi = T(reference)
		psi.increment(reference)
		E = dot(psi,H(psi)) / dot(psi,psi)
		print("Evar = {}  (error = {})".format(E,E+2.8717891505016))
	#
	phase = +1
	for i in range(3):
		phase *= -1
		accumulate.increment(intermediate,phase)
		if i==2:  break
		intermediate = T(intermediate)
	accumulate = R(accumulate)
	dT = SD_excitations(reference, accumulate, cisd_strings, occ_orbs, vrt_orbs)
	T.increment(dT)
