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
#
# ANALYTICAL SOLUTION = 2.398435656980412
# CISD Results:
#              5 x 5 => 2.3984416499
#              7 x 7 => 2.39843571482
#

import sys
import time
import copy
import numpy as np
from qode.util import parallel, output, textlog, indent

# Config and state infrastructure, SCF and CCSD modules
from qode.fermion_field import state
from qode.SCF_diag_full import hf as scf_routine
from qode.many_body     import CCSD
from Applications.component_tests.oscillator_ccsd.old_qode.generate_ho_config import diag_rec_h0_mat, diag_rec_v_matrix, generate_main, get_lowest_analytical, get_rec_obj

# Oscillator System Imports
from qode.coupled_oscillators            import oscillator_system as osys
from qode.coupled_oscillators            import hamiltonian_class
from qode.coupled_oscillators.analytical import analytical_solution

# H and T operators, and their manipulations
from hamiltonian import hamiltonian, inv_fock
from t_operator  import T_operator, T_space_traits
from bch         import BakerCampbellHausdorff

def get_analytical_solution( mol_obj ):
    required_states = 5
    k_eff, states_configure, E_analytical, vectors = \
        analytical_solution.analytical_main(
                                            mol_obj.list_of_masses() ,
                                            mol_obj.get_k_mat()      ,
                                            required_states,
                                            1.0,
                                            1.0,
                                            1.0
                                            )
    #print("Lowest Analytical Solution = ", sorted( E_analytical )[0] )
    return sorted( E_analytical )[0]

def main(h_mat, V_mat, rec_num_states, textlog, resources):
	out = output(log=textlog)

	out.log("Parsing input ...")
	cisd_strings = generate_main(rec_num_states)
	out.log(indent("Recoupling Number of States =", rec_num_states))

	out.log("Building H ...")
	H = hamiltonian(h_mat, V_mat, rec_num_states)

	out.log("Getting orbital energy preconditioner ...")
	invF = inv_fock(h_mat, rec_num_states)

	out.log("Initializing T ...")
	t_amp_obj =  state.state(cisd_strings)	# Null vector
	t_amp_obj.update_coeffs( np.zeros(( len(cisd_strings) )) )
	T = T_operator(t_amp_obj, rec_num_states)

	BCH = BakerCampbellHausdorff(cisd_strings, resources)
	out.log("Baker-Campbell-Hausdorff module initialized.")

	out.log("Perform coupled-cluster iterations.")
	t_ccsd_start = time.time()
	E = CCSD(H, T, invF, BCH, T_space_traits, out.log.sublog())
	t_ccsd_end = time.time()

	out.log("CCSD Energy =", E, 'Hartree')
	out.log("CCSD TIME {}".format(t_ccsd_end - t_ccsd_start))

	return out(energy=E)



if __name__ == "__main__":
	out = output(log=textlog(echo=True))
	resources = parallel.resources(int(sys.argv[1]))
	# 
	rec_num_states = [10, 10]
	#
	ho1 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho2 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	#
	ho3 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho4 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	#
	coupling_mat = [[0.0, 0.3],[0.3,0.0]]
	#
	mol1 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho1,ho2] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	#
	mol2 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho3,ho4] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	#
	ext_k_mat = [[0.0,0.1],[0.1,0.0]]
	mol_obj = osys.OscillatorSystem( osys.group_of_oscillators( recursive_list =  [mol1, mol2],
                                                                 coupling     =  ext_k_mat ) )
	#
	rec_obj = get_rec_obj( mol_obj, rec_num_states )
	eig1, eig2 = rec_obj.get_sorted_eigen_energy()
	v_mat_list = rec_obj.get_pair_dipole_mat()
	#
	# Building h_mat and V_mat objects
	h_mat = diag_rec_h0_mat([eig1, eig2])
	V_mat = diag_rec_v_matrix(rec_num_states, v_mat_list)
	#
	out.log('### COUPLED-CLUSTER CALCULATION ###\n')
	#
	t1 = time.time()
	#
	out.ccsd = \
		main(h_mat, V_mat, rec_num_states, textlog=out.log.sublog(), resources=resources)
	#
	t2 = time.time()
	#
	out.log('\n### RESULTS ###\n')
	out.log("TOTAL ENERGY =", out.ccsd.energy)
	out.log("TOTAL TIME   =", t2 - t1, "seconds\n")
	out.log.write('log.out')
	out.log("------------------------------------")
	out.log("ANALYTICAL SOLUTION =",get_lowest_analytical(mol_obj))
	out.log("------------------------------------")
