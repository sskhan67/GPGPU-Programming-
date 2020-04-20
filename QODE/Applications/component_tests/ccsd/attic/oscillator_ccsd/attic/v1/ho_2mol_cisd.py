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
# System Imports
import copy
import numpy as np
from multiprocessing import Pool

# Oscillator System Imports
from qode.coupled_oscillators import oscillator_system as osys
from qode.coupled_oscillators import hamiltonian_class
from qode.coupled_oscillators.analytical import analytical_solution

# config/state Imports
from qode.fermion_field import state, config


# Generation of CISD Vector Imports
from generate_ho_config import generate_main, get_lowest_analytical, get_rec_obj, \
                                                                  diag_rec_h0_mat, diag_rec_v_matrix

# H Operator Imports
from osc_hamiltonian_operator import h_wrapper



def two_mol_CISD_main(state_obj, H, V, rec_num_states):
	# state_obj is a state.state class object with CISD configs 
	# H, V are two Matrix Classes holding fragment eigenvalues and coupling matrix values
	# rec_num_states is a list of max number of states used in this calculation for each fragment.
	#
	#
	# Make Parallel Input List,
	input_list = []
	occ_string = state_obj.value()[0]
	cisd_mat_dim = len(occ_string)
	for j in range(cisd_mat_dim):
		state_j = state.state(occ_string)
		new_coeffs = np.zeros((cisd_mat_dim))
		new_coeffs[j] = 1.0
		state_j.update_coeffs( new_coeffs )
		input_list += [ [state_j, H, V, rec_num_states] ]
	#
	# Send them to map function,
	myPool = Pool(25)
	generated_states = myPool.map(h_wrapper, input_list)
	myPool.terminate()
	#
	# Loop Over CISD states to get CISD Matrix,
	H_2nd = []
	#
	#
	for i in range(cisd_mat_dim):   # First Dimension
		state_i = state.state(occ_string)
		new_coeffs = np.zeros((cisd_mat_dim))
		new_coeffs[i] = 1.0
		state_i.update_coeffs(new_coeffs)
		for j in range(cisd_mat_dim):   # Second Dimension
			H_2nd += [ state.dot( state_i, generated_states[j] ) ]
			# print("< %d | H | %d > = %f" %( i,j, H_2nd[-1] ) )
	#
	# Reshape into Matrix and do diagonaliztion,
	H_2nd = np.array(H_2nd)
	H_2nd = np.reshape(H_2nd, (cisd_mat_dim, cisd_mat_dim))
	H_2nd = np.matrix(H_2nd)
	E,U = np.linalg.eigh(H_2nd)
	print("ANALYTICAL SOLUTION =",get_lowest_analytical(mol_obj))
	print("Lowest Eigen Value  =", E[0])
	# 
	# Returns Lowest Eigenvalue.
	return E[0]









if __name__ == "__main__":
	#
	# CISD Results:
	# ANALYTICAL SOLUTION = 2.398435656980412
	#              5 x 5 => 2.3984416499
	#              7 x 7 => 2.39843571482
	# 10x10 =>
	rec_num_states = [5,5]
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
	#
	ext_k_mat = [[0.0,0.1],[0.1,0.0]]
	mol_obj = osys.OscillatorSystem( osys.group_of_oscillators( recursive_list =  [mol1, mol2],
                                                                 coupling      =  ext_k_mat ) )
	#
	#
	#
	rec_obj = get_rec_obj( mol_obj, rec_num_states )
	eig1, eig2 = rec_obj.get_sorted_eigen_energy()
	#
	v_mat_list = rec_obj.get_pair_dipole_mat()
	#
	#
	# H here is the h matrix object
	H = diag_rec_h0_mat([eig1, eig2])
	#
	#
	# V here is the coupling matrix object
	V = diag_rec_v_matrix(rec_num_states, v_mat_list)
	#
	#
	#
	ho_configs = generate_main(rec_num_states)
	state_obj  = state.state(ho_configs)
	# NOW COMES THE CISD FUNCTION CALL.
	#
	#
	E_CISD = two_mol_CISD_main(state_obj, H, V, rec_num_states)
	#




