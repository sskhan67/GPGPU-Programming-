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
import sys
import time
import signal
import traceback
import copy
import numpy
import ctypes as ct
from qode.util import parallel, output, textlog, indent

# CCSD module
import ccsd

# Oscillator System Imports
from qode.coupled_oscillators.CI_util.generate_ho_config import get_lowest_analytical, get_rec_obj, diag_rec_h0_mat, diag_rec_v_matrix
from qode.coupled_oscillators            import oscillator_system as osys
from qode.coupled_oscillators            import hamiltonian_class
from qode.coupled_oscillators.analytical import analytical_solution


def timeout_handler(signum, frame):
    print("Single Job Max Time Exceeded!")
    raise signal.ItimerError("Single_Job_Time_Out")


def get_primitive_mol_obj(
                num_oscillator,
                num_molecule,
                positions
                ):

	# k list for one molecule
	k_list = [ round(  i/(num_oscillator-1)+1.0, 4 ) for i in range(num_oscillator) ]
	    
	# internal couplings for one molecule
	coupling_mat = [[ 0.0 for i in range(num_oscillator) ] for j in range(num_oscillator) ]

	# fill in internal couplings
	for i in range( num_oscillator ):
		for j in range( i+1, num_oscillator ):
			coupling_mat[i][j] = round( ( k_list[j] - k_list[i] )/ 3.0 , 4 )

	# external couplings for "num_molecule" molecules.
	ext_k_mat = [[ 0.0 for i in range(num_molecule * num_oscillator) ] for j in range(num_molecule * num_oscillator) ]

	for i in range(num_molecule):
		for x in range(num_oscillator):
			for y in range(x+1, num_oscillator):
				ext_k_mat[i*num_oscillator+x][i*num_oscillator+y] = coupling_mat[x][y]

	for i in range(num_molecule):
		for j in range(num_molecule):
			if i != j:
				for x in range(num_oscillator):
					for y in range(num_oscillator):
						ext_k_mat[i*num_oscillator + x][j*num_oscillator +y] = -2.0/ abs(positions[i] - positions[j])**3
	# for line in coupling_mat:
	# 	print(line)      
	# for line in ext_k_mat:
	# 	for item in line:
	# 		print("%.3e" %(item), end='\t')
	# 	print()

	oscillators = [ osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = k_list[i] ) )\
	                      for i in range(num_oscillator) ]

	fragments = []
	for i in range(num_molecule):
		for j in range(num_oscillator):
			fragments += [ osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ oscillators[j] ] ),
	                                                                        coupling       = [1.0]  ) ) ]

	        
	molecule = osys.OscillatorSystem( osys.group_of_oscillators( recursive_list = fragments,
	                                                             coupling       = ext_k_mat ) )
	        
	 
	# print(molecule.get_k_mat())
	# print(molecule.top_level_groups())
	# print(molecule.top_level_coupling_mat())
	return molecule


def main(n_ho, n_state, n_mol, mol_obj, analytical_value):
	print("+++++++++++++++++++++ %d ho %d mol %d states  ++++++++++++++++++++++++++++++++++++++" %(n_ho, n_mol, n_state))
	out = output(log=textlog(echo=True))
	resources = parallel.resources(int(sys.argv[1]))
	 
	t1 = time.time()

	out.log('\n### SET UP OSCILLATOR SYSTEMS ###\n')

	rec_num_states = [n_state for i in range(n_mol*n_ho)]
	rec_obj = get_rec_obj( mol_obj, rec_num_states )
	
	all_mol_eigenvalues = rec_obj.get_sorted_eigen_energy()
	v_mat_list = rec_obj.get_pair_dipole_mat()
	

	# raise RuntimeError
	out.log("Parsing CC input ...")
	# out.log(indent("total number of orbitals = ", n_orbs))
	# out.log(indent("occupied orbital indices = ", occ_orbs))
	# out.log(indent("virtual  orbital indices = ", vrt_orbs))

	# Building h_mat and V_mat objects
	out.log('### COUPLED-CLUSTER CALCULATION ###\n')

	# out.log(indent("occupied orbital indices = ", occ_orbs))
	# out.log(indent("virtual  orbital indices = ", vrt_orbs))


	F_mat = all_mol_eigenvalues
	V_mat = v_mat_list
	# out.ccsd = ccsd.main( [], V_mat, F_mat, rec_num_states, textlog=out.log.sublog(), resources=resources )
	max_time = 3600*3 
	signal.signal(signal.SIGALRM, timeout_handler)
	signal.alarm(max_time)

	try:
	    out.ccsd = ccsd.main( [], V_mat, F_mat, rec_num_states, textlog=out.log.sublog(), resources=resources )
	except signal.ItimerError:
	    print("Calculation TimeOut!")
	    exc_type, exc_value, exc_traceback = sys.exc_info()
	    print(traceback.format_exception(exc_type, exc_value,exc_traceback))
	    return False
	except:
	    print("UNKNOWN ERROR!")
	    exc_type, exc_value, exc_traceback = sys.exc_info()
	    print(traceback.format_exception(exc_type, exc_value,exc_traceback))
	    return False

	t2 = time.time()

	# analytical_value = get_lowest_analytical(mol_obj)
	out.log('\n### RESULTS ###\n')
	out.log("TOTAL ENERGY =", out.ccsd.energy)
	out.log("TOTAL TIME   =", t2 - t1, "seconds\n")
	out.log("------------------------------------")
	out.log("ANALYTICAL SOLUTION =", analytical_value)
	out.log("------------------------------------")
	out.log.write('log.out')
	return True

def analytical_with_timer(mol_obj):
	max_time = 3600*3 
	signal.signal(signal.SIGALRM, timeout_handler)
	signal.alarm(max_time)
	try:
		return get_lowest_analytical(mol_obj)
	except signal.ItimerError:
	    print("Calculation TimeOut!")
	    exc_type, exc_value, exc_traceback = sys.exc_info()
	    print(traceback.format_exception(exc_type, exc_value,exc_traceback))
	    return None
	except:
	    print("UNKNOWN ERROR!")
	    exc_type, exc_value, exc_traceback = sys.exc_info()
	    print(traceback.format_exception(exc_type, exc_value,exc_traceback))
	    return None





if __name__ == "__main__":
	import traceback
	n_ho = 8
	for n_mol, n_state_max in [ [2,10],[3,10],[4,10],[5,10],[6,10],[7,10],[8,10],[9,10],[10,9],
	                            [11,9],[12,9],[13,8],[14,8],[15,7], 
	                            [16,7],
                                [17,7],
                                [18,7],[19,6],
	                            [20,5],[21,5],[22,5],[23,5],
                                [24,5],
		                        [25,3], [26,3], [27,3],[28,3],[29,3],
                                [30,3]]:
		positions = [ 5.0 * i for i in range(n_mol)]
		t_start_mol = time.time()
		mol_obj = get_primitive_mol_obj(
				n_ho,
                n_mol,
                positions
                )
		t_end_mol = time.time()
		print("Building Molecular Object Time = %f sec" %(t_end_mol- t_start_mol))

		t_start_ana = time.time()
		analytical_value = get_lowest_analytical(mol_obj)
		t_end_ana = time.time()

		# max_time = 3600*3 
		# signal.signal(signal.SIGALRM, timeout_handler)
		# signal.alarm(max_time)

		if analytical_value != None:
			print("Analytical Solution      = %f" %(analytical_value))
			print("Analytical Solution Time = %f sec" %(t_end_ana- t_start_ana))
			n_state = 3
			while n_state <= n_state_max and main(n_ho, n_state, n_mol, mol_obj, analytical_value):
				# print("INCREASING NUM STATES TO", n_state + 1)
				n_state += 1
			# Some Output to Help to understand...
			if n_state < n_state_max:
				n_state += 1
				while n_state <= n_state_max:
					print("+++++++++++++++++++++ %d ho %d mol %d states  ++++++++++++++++++++++++++++++++++++++\nFAILED!!" %(n_ho, n_mol, n_state))
					n_state += 1

		else:
			print("+++++++++++++++++++++ %d ho %d mol           ++++++++++++++++++++++++++++++++++++++\nFAILED!!" %(n_ho, n_mol))





