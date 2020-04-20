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
from qode.util import parallel, output, indent, textlog

# operators CCSD operators and modules
from qode.many_body import CCSD
from qode.many_body.hierarchical_fluctuations.excitonic import H_operator, T_operator, space_traits, BakerCampbellHausdorff, inv_frag_eigvals

# Oscillator System Imports
from qode.coupled_oscillators                            import oscillator_system as osys
from qode.coupled_oscillators.CI_util.generate_ho_config import get_lowest_analytical, get_rec_obj

# Hamiltonian
import hamiltonian



def get_mol_obj(num_oscillator, num_molecule, positions, INTER_ON=True):
    # k list for one molecule
    k_list = [ round(1+i/(num_oscillator-1.),4) for i in range(num_oscillator) ]
    # internal couplings for one molecule
    coupling_mat = [[ 0. for _1 in range(num_oscillator) ] for _2 in range(num_oscillator) ]
    for i in range(num_oscillator):
        for j in range(i+1,num_oscillator):  coupling_mat[i][j] = round( (k_list[j]-k_list[i])/3, 4 )
    # external couplings for "num_molecule" molecules.
    ext_k_mat = [[ 0. for _1 in range(num_molecule) ] for _2 in range(num_molecule) ]
    if INTER_ON:
	    for i in range(num_molecule):
	        for j in range(num_molecule):
	            if i != j:  ext_k_mat[i][j] = -2.0 / abs(positions[i]-positions[j])**3
    # final setup
    oscillators = [ osys.OscillatorSystem( osys.primitive_oscillator(m=1.,k=k) ) for k in k_list ]
    fragments   = [ osys.OscillatorSystem( osys.group_of_oscillators(recursive_list=copy.deepcopy(oscillators),coupling=copy.deepcopy(coupling_mat)) ) for _ in range(num_molecule) ]
    molecule    =   osys.OscillatorSystem( osys.group_of_oscillators(recursive_list=fragments,coupling=ext_k_mat) )
    return molecule

def build(separation, n_mol, n_ho):
	positions = [ separation*i for i in range(n_mol)]
	mol_obj        = get_mol_obj(n_ho, n_mol, positions, INTER_ON=True)
	mol_obj_no_int = get_mol_obj(n_ho, n_mol, positions, INTER_ON=False)
	return mol_obj, mol_obj_no_int

def do_analytical(mol_obj, mol_obj_no_int, out):
	analytical_value        = get_lowest_analytical(mol_obj)
	analytical_value_no_int = get_lowest_analytical(mol_obj_no_int)
	analytical_int_E = abs(analytical_value - analytical_value_no_int)
	out.log("Analytical Energy (with interaction) =", analytical_value)
	out.log("Analytical Energy (NO   interaction) =", analytical_value_no_int)
	out.log("Analytical Interaction Energy        =", analytical_int_E)
	return analytical_value



def main(mol_obj, rec_num_states, out, analytical_value, resources):
	t1 = time.time()

	out.log("+++++++++++++++++++++ %d ho %d mol %d states  ++++++++++++++++++++++++++++++++++++++" %(n_ho, n_mol, n_state))
	out.log('\n### SET UP OSCILLATOR SYSTEMS ###\n')
	rec_obj = get_rec_obj( mol_obj, rec_num_states )
	F_mat = rec_obj.get_sorted_eigen_energy()
	V_mat = rec_obj.get_pair_dipole_mat()
	out.log("Parsing CC input ...")
	out.log('### COUPLED-CLUSTER CALCULATION ###\n')

	out.ccsd = output(log=out.log.sublog())

	out.ccsd.log("Building blocked H ...")
	H = H_operator(rec_num_states)
	hamiltonian.load(H, rec_num_states, F_mat, V_mat)

	out.ccsd.log("Getting orbital energy preconditioner ...")
	invF = inv_frag_eigvals(F_mat, rec_num_states)

	out.ccsd.log("Initializing T ...")
	T = T_operator(rec_num_states)		# Initialized to zero

	BCH = BakerCampbellHausdorff( rec_num_states, resources )
	out.ccsd.log("Baker-Campbell-Hausdorff module initialized.")

	out.ccsd.log("Perform coupled-cluster iterations.")
	t_ccsd_start = time.time()
	E0, H.K[0] = H.K[0], 0.
	out.ccsd.energy = E0 + CCSD(H, T, invF, BCH, space_traits, out.ccsd.log.sublog(), thresh=1e-15)
	t_ccsd_end = time.time()

	out.ccsd.log("CCSD Energy =", out.ccsd.energy, 'Hartree')
	out.ccsd.log("CCSD TIME {}".format(t_ccsd_end - t_ccsd_start))

	t2 = time.time()

	out.log('\n### RESULTS ###\n')
	out.log("TOTAL ENERGY =", out.ccsd.energy)
	out.log("TOTAL TIME   =", t2 - t1, "seconds\n")
	out.log("------------------------------------")
	out.log("ANALYTICAL SOLUTION =", analytical_value)
	out.log("------------------------------------")





if __name__ == "__main__":
	resources, out  = parallel.resources(int(sys.argv[2])), output(log=textlog(echo=True))
	separation, n_ho = float(sys.argv[1]), 8
	for n_mol,max_states in [         [2,20], [3,20], [4,20], [5,20], [6,20], [7,20], [8,20], [9,20],[10,20],
	                         [11,10],[12,10],[13,10],[14,10],[15,10],[16,10],[17,10],[18,10],[19,10],[20,10],
	                         [21,10],[22,10],[23,10],[24,10],[25,10],[26,10],[27,10],[28,10],[29,10],[30,10]]:
		mol_obj, mol_obj_no_int = build(separation, n_mol, n_ho)
		analytical_value = do_analytical(mol_obj, mol_obj_no_int, out)
		for n_state in range(3,max_states+1):
			main(mol_obj, [n_state for i in range(n_mol)], out, analytical_value, resources)
	out.log.write('excitonic_{:.1f}.out'.format(separation))










