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
# utilities
import time
from qode.util import output, indent

# operators, SCF and CCSD modules
from qode                               import fermion_field
from qode.many_body.hierarchical_fluctuations.local_operator      import space_traits, BakerCampbellHausdorff
from qode.many_body                     import CCSD
from qode.many_body.hierarchical_fluctuations.local_operator.local_operator import T_operator

# Hamiltonian
from hamiltonian import load_values_into_blocks
from qode.many_body.hierarchical_fluctuations.local_operator.space_traits import orbital_energy



def main(h_mat, V_mat, F_mat, rec_num_states, textlog, resources):
	out = output(log=textlog)

	# out.log("Partioning fluctuations (excitations, flat, deexcitations) ...")

	out.log("Building blocked H ...")
	H = load_values_into_blocks(rec_num_states, F_mat, V_mat)
	E0 = H.K[0]
	H.K[0] = 0.0
	# print("H ------------------------------------------------------")
	# H_dict = H.__dict__
	# for key,values in H_dict.items():
	# 	if '_' not in key:
	# 		if values != None:
	# 			print(key)
	# 			print(values)
	# print("H ------------------------------------------------------")

	out.log("Getting orbital energy preconditioner ...")
	invF = orbital_energy(F_mat, rec_num_states)

	out.log("Initializing T ...")
	T = T_operator(rec_num_states)		# Initialized to zero

	BCH = BakerCampbellHausdorff( rec_num_states, resources )
	out.log("Baker-Campbell-Hausdorff module initialized.")

	out.log("Perform coupled-cluster iterations.")
	t_ccsd_start = time.time()
	E = CCSD(H, T, invF, BCH, space_traits, out.log.sublog())
	t_ccsd_end = time.time()

	E += E0
	out.log("CCSD Energy =", E, 'Hartree')
	out.log("CCSD TIME {}".format(t_ccsd_end - t_ccsd_start))

	return out(energy=E)
