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
import copy
import time
import numpy as np
from math import sqrt
from multiprocessing import Pool

# Config/state imports
from qode.fermion_field import occ_strings, state

# Hamiltonian imports
from Applications.ccsd.elec_ccsd.elec_hamiltonian import h_wrapper




def cisd_main(h_mat, V_mat, num_alpha_elec, num_beta_elec, num_spin_orb):
	#
	#########
	occ_orbs = []
	vrt_orbs = []
	for p in range(num_alpha_elec):  occ_orbs += [p]
	for p in range(num_beta_elec):   occ_orbs += [p+(num_spin_orb//2)]
	for p in range(num_spin_orb):
		if p not in occ_orbs:  vrt_orbs += [p]
	#########
	cisd_strings = occ_strings.CISD(occ_orbs,vrt_orbs)
	state_obj    = state.state(cisd_strings)
	#
	#
	#  Now Do CISD Code.
	#
	num_config = len(cisd_strings)
	#
	#
	input_list = []
	for i in range(num_config):
		state_obj_i = state.state(cisd_strings,i)
		input_list += [[state_obj_i, h_mat, V_mat]]
	#
	myPool = Pool(25)
	generated_states = myPool.map(h_wrapper, input_list)
	myPool.terminate()
	#
	#
	#
	H_2nd = []
	#
	#
	#
	for i in range(num_config):
		state_i = state.state(cisd_strings,i)
		#
		#
		for j in range(num_config):
			H_2nd += [ state.dot( state_i, generated_states[j] ) ]
			print("< %d | H | %d > = %f" %( i,j, H_2nd[-1] ) )
	#
	#
	print(H_2nd)
	H_2nd = np.array(H_2nd)
	H_2nd = np.reshape(H_2nd, (num_config,num_config))
	H_2nd = np.matrix(H_2nd)
	#
	#
	# np.save('SecondQuantizedExplicitState_H_mat.npy', H_2nd)
	#
	#
	E,U = np.linalg.eigh(H_2nd)
	print("Lowest Eigen Value =", E[0])
	return E[0]





if __name__ == "__main__":
	import qode.SCF_diag_full.hf as scf_routine
	#
	# Helium Atom:  CCSD/6-31G uncontracted
	# SCF energy in the final basis set = -2.8551604262  Hartree (Qchem4.4)
	# CCSD total energy                 = -2.87178917    Hartree (Qchem4.4)
	# 
	#
	inputs = ['he']
	for item in inputs:
		t1 = time.time()
		input_file = item + '.in'
		E_0, KE, NAE, NRE, ERE, F_mat, T_mat, N_mat, V_mat, vHF, num_alpha_elec, num_beta_elec, num_spin_orb = scf_routine.main(input_file, 1e-10)
		h_mat = T_mat + N_mat
		E_CISD = cisd_main(h_mat, V_mat, num_alpha_elec, num_beta_elec, num_spin_orb)
		# f = open(item + '.out', 'w')
		# f.write('E_CCSD = ' + str(E_CCSD))
		# f.close()
		t2 = time.time()
		print("CISD TOTAL TIME =", t2 - t1, "seconds")
