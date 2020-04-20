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
import copy
import numpy
import ctypes as ct
from qode.util import parallel, output, textlog, indent

# CCSD module
from local import ccsd

# Oscillator System Imports
from attic.oscillator_ccsd.old_qode.generate_ho_config import get_lowest_analytical, get_rec_obj, diag_rec_h0_mat, diag_rec_v_matrix
from qode.coupled_oscillators            import oscillator_system as osys
from qode.coupled_oscillators            import hamiltonian_class
from qode.coupled_oscillators.analytical import analytical_solution




def get_mol_obj(
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
            coupling_mat[i][j] = round( ( k_list[i] - k_list[j] )/ 3.0 , 4 )

    # external couplings for "num_molecule" molecules.
    ext_k_mat = [[ 0.0 for i in range(num_molecule) ] for j in range(num_molecule) ]
    
    for i in range(num_molecule):
        for j in range(num_molecule):
            if i != j:
                ext_k_mat[i][j] = -2.0/ abs(positions[i] - positions[j])**3

                    
    oscillators = [ osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = k_list[i] ) )\
                          for i in range(num_oscillator) ]

    fragments = [ osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( oscillators  ),
                                                                   coupling       = copy.deepcopy( coupling_mat )  ) )
                          for i in range(num_molecule) ]
            
    molecule = osys.OscillatorSystem( osys.group_of_oscillators( recursive_list = fragments,
                                                                 coupling       = ext_k_mat ) )
            
     
    # print(molecule.get_k_mat())
    # print(molecule.top_level_groups())
    # print(molecule.top_level_coupling_mat())
    return molecule


def main(n_ho, n_state, n_mol):
	out = output(log=textlog(echo=True))
	resources = parallel.resources(int(sys.argv[1]))
	 
	t1 = time.time()

	out.log('\n### SET UP OSCILLATOR SYSTEMS ###\n')


	rec_num_states = [n_state for i in range(n_mol)]
	positions = [ 10.0 * i for i in range(n_mol)]


	mol_obj = get_mol_obj(
				n_ho,
                n_mol,
                positions
                )


	rec_obj = get_rec_obj( mol_obj, rec_num_states )
	
	all_mol_eigenvalues = rec_obj.get_sorted_eigen_energy()
	v_mat_list = rec_obj.get_pair_dipole_mat()
	
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
	out.ccsd = ccsd.main( [], V_mat, F_mat, rec_num_states, textlog=out.log.sublog(), resources=resources )

	t2 = time.time()
	out.log('\n### RESULTS ###\n')
	out.log("TOTAL ENERGY =", out.ccsd.energy)
	out.log("TOTAL TIME   =", t2 - t1, "seconds\n")
	out.log("------------------------------------")
	out.log("ANALYTICAL SOLUTION =",get_lowest_analytical(mol_obj))
	out.log("------------------------------------")
	out.log.write('log.out') 


if __name__ == "__main__":
	import traceback
	n_ho = 8
	for n_mol in [5]: #[10,35,40,45,50,55,60,65,70,75,80,85,90,95,100]:
		for n_state in [10]:
		# for n_mol, n_state in [[3,5]]:
			try:
				print("+++++++++++++++++++++ %d ho %d mol %d states  ++++++++++++++++++++++++++++++++++++++" %(n_ho, n_mol, n_state))
				main(n_ho, n_state, n_mol)
				print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
			except:
				print("----------------------------------------------------------------------------")
				print("FAILED!\n----------------------------------------------------------------------------")
				exc_type, exc_value, exc_traceback = sys.exc_info()
				traceback.print_tb(exc_traceback)
				print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")




