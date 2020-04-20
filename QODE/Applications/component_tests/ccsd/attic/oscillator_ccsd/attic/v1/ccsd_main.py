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

import time
import copy
from qode import util
import numpy as np

# Oscillator System Imports
from qode.coupled_oscillators import oscillator_system as osys
from qode.coupled_oscillators import hamiltonian_class
from qode.coupled_oscillators.analytical import analytical_solution

# Config and state infrastructure
from qode.fermion_field import state, config
from generate_ho_config import diag_rec_h0_mat, diag_rec_v_matrix, generate_main, get_lowest_analytical, get_rec_obj

# CCSD modules
from qode.many_body import CCSD

# H and T operators
#from Applications.ccsd.elec_ccsd.hamiltonian_function import hamiltonian, inv_fock
from osc_hamiltonian_operator   import inv_fock
from osc_hamiltonian_operator   import hamiltonian
from osc_t_operator             import T_operator




# def main(h_mat, V_mat, rec_num_states, num_max_subspace=5, num_itr_start_diis=1):

def main(h_mat, V_mat, rec_num_states, num_max_subspace, num_itr_start_diis, textlog):
	output = util.output()
	output.log = textlog
	output.log_indent = util.indented(textlog)

	output.log("Parsing input ...")
	# occ_orbs = []
	# vrt_orbs = []
	# for p in range(num_alpha_elec):  occ_orbs += [p]
	# for p in range(num_beta_elec):   occ_orbs += [p+(num_spin_orb//2)]
	# for p in range(num_spin_orb):
	# 	if p not in occ_orbs:  vrt_orbs += [p]
	cisd_strings = generate_main(rec_num_states)

	# output.log_indent("occupied orbital indices = ", occ_orbs)
	# output.log_indent("virtual  orbital indices = ", vrt_orbs)
	output.log_indent("Recoupling Number of States =", rec_num_states)




	output.log("Building H ...")
	#H = hamiltonian(h_mat, V_mat, cisd_strings)								# hamiltonian_simple module
	# H = hamiltonian( h_mat, V_mat, cisd_strings, list(reversed(occ_orbs)), list(reversed(vrt_orbs)) )	# hamiltonian_normal module
	H = hamiltonian(h_mat, V_mat, rec_num_states)


	output.log("Getting orbital energy preconditioner ...")
	invF = inv_fock(h_mat, rec_num_states)
	

	output.log("Initializing T ...")
	t_amp_obj =  state.state(cisd_strings)
	t_amp_obj.update_coeffs( np.zeros(( len(cisd_strings) )) )
	# offsets = ( 1, len(occ_strings.CIS(occ_orbs,vrt_orbs)) )
	# T = T_operator(t_amp_obj, occ_orbs, vrt_orbs, offsets)
	T = T_operator(t_amp_obj, rec_num_states)



	output.log("Begin coupled-cluster iterations.")
	t_ccsd_start = time.time()
	E_nuclear_repul = 0.0  # No Nuclear Repulsion in this case
	output.ccsd_solver = \
		CCSD(cisd_strings, H, T, E_nuclear_repul, invF, num_max_subspace, num_itr_start_diis, textlog=output.log.sublog()).iterate_to_convergence()
	t_ccsd_end = time.time()

	output.log("CCSD converged at Cycle", output.ccsd_solver.cycle)
	output.log("CCSD Energy =", output.ccsd_solver.E, 'Hartree')
	output.log("CCSD TIME {}".format(t_ccsd_end - t_ccsd_start))

	return output(energy=output.ccsd_solver.E)





if __name__ == "__main__":
	output = util.output()
	output.log = util.textlog(echo=True)
	# 
	# CISD Results:
	# ANALYTICAL SOLUTION = 2.398435656980412
	#              5 x 5 => 2.3984416499
	#              7 x 7 => 2.39843571482
	rec_num_states = [5,5]
	#
	#
	#
	ho1 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho2 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	#
	#
	ho3 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho4 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	#
	#
	coupling_mat = [[0.0, 0.3],[0.3,0.0]]
	#
	#
	mol1 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho1,ho2] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	#
	#
	mol2 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho3,ho4] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	#
	#
	ext_k_mat = [[0.0,0.1],[0.1,0.0]]
	mol_obj = osys.OscillatorSystem( osys.group_of_oscillators( recursive_list =  [mol1, mol2],
                                                                 coupling     =  ext_k_mat ) )
	#
	#
	rec_obj = get_rec_obj( mol_obj, rec_num_states )
	eig1, eig2 = rec_obj.get_sorted_eigen_energy()
	v_mat_list = rec_obj.get_pair_dipole_mat()
	#
	# Building h_mat and V_mat objects
	h_mat = diag_rec_h0_mat([eig1, eig2])
	V_mat = diag_rec_v_matrix(rec_num_states, v_mat_list)
	#
	#
	output.log('### COUPLED-CLUSTER CALCULATION ###\n')
	#
	#
	t1 = time.time()

	output = main(h_mat, V_mat, rec_num_states, num_max_subspace=5, num_itr_start_diis=1, textlog=output.log.sublog())
	
	t2 = time.time()

	output.log('\n### RESULTS ###\n')
	output.log("TOTAL ENERGY =", output.energy)
	output.log("TOTAL TIME   =", t2 - t1, "seconds")
	# output.log.write(item + '.out')
	#
	# np.savetxt("HO_CCSD" + time.ctime() + ".txt", np.matrix([[E_CCSD]]))
	#
	output.log("------------------------------------")
	output.log("ANALYTICAL SOLUTION =",get_lowest_analytical(mol_obj))
	output.log("------------------------------------")
	#
