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
import copy
from qode.util import parallel, output, indent, textlog

# operators CCSD operators and modules
from qode.many_body import CCSD
#from qode.many_body.hierarchical_fluctuations.excitonic import H_operator, T_operator, space_traits, BakerCampbellHausdorff, inv_frag_eigvals
from qode.many_body.hierarchical_fluctuations.local_operator import H_operator, T_operator, space_traits, BakerCampbellHausdorff
from qode.many_body.hierarchical_fluctuations.local_operator import orbital_energy as inv_frag_eigvals

# Oscillator System Imports
from qode.coupled_oscillators                            import oscillator_system as osys
from qode.coupled_oscillators.CI_util.generate_ho_config import get_lowest_analytical, get_rec_obj

# Hamiltonian
import hamiltonian

# maximum length of a job in seconds (3 hours)
def timeout_handler(signum, frame):   raise signal.ItimerError()
signal.signal(signal.SIGALRM, timeout_handler)
max_time = 3*3600



def get_primitive_mol_obj(num_oscillator, num_molecule, positions):
	# k list for one molecule
	k_list = [ round(1+i/(num_oscillator-1.),4) for i in range(num_oscillator) ]
	# internal couplings for one molecule
	coupling_mat = [[ 0. for _1 in range(num_oscillator) ] for _2 in range(num_oscillator) ]
	for i in range(num_oscillator):
		for j in range(i+1,num_oscillator):  coupling_mat[i][j] = round( (k_list[j]-k_list[i])/3, 4 )
	# external couplings for "num_molecule" molecules.
	ext_k_mat = [[ 0. for _1 in range(num_molecule*num_oscillator) ] for _2 in range(num_molecule*num_oscillator) ]
	for i in range(num_molecule):
		for x in range(num_oscillator):
			for y in range(x+1, num_oscillator):
				ext_k_mat[i*num_oscillator+x][i*num_oscillator+y] = coupling_mat[x][y]
	for i in range(num_molecule):
		for j in range(num_molecule):
			if i != j:
				for x in range(num_oscillator):
					for y in range(num_oscillator):
						ext_k_mat[i*num_oscillator+x][j*num_oscillator+y] = -2.0 / abs(positions[i]-positions[j])**3
	# final setup
	oscillators = [ osys.OscillatorSystem( osys.primitive_oscillator(m=1.,k=k) ) for k in k_list ]
	fragments = []
	for _ in range(num_molecule):
		fragments += [ osys.OscillatorSystem( osys.group_of_oscillators(recursive_list=copy.deepcopy([osc]), coupling=[1.0])) for osc in oscillators ]
	molecule    =   osys.OscillatorSystem( osys.group_of_oscillators(recursive_list=fragments,coupling=ext_k_mat) )
	return molecule

def build(separation, n_mol, n_ho):
	positions = [ separation*i for i in range(n_mol)]
	t_start_mol = time.time()
	mol_obj = get_primitive_mol_obj(n_ho, n_mol, positions)
	t_end_mol = time.time()
	out.log("Building Molecular Object Time = %f sec" %(t_end_mol- t_start_mol))
	return mol_obj

def do_analytical(mol_obj, out):
	t_start_ana = time.time()
	analytical_value = get_lowest_analytical(mol_obj)
	t_end_ana = time.time()
	if analytical_value != None:
		out.log("Analytical Solution      = %f" %(analytical_value))
		out.log("Analytical Solution Time = %f sec" %(t_end_ana- t_start_ana))
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
	signal.alarm(max_time)
	try:
		E0, H.K[0] = H.K[0], 0.
		out.ccsd.energy = E0 + CCSD(H, T, invF, BCH, space_traits, out.ccsd.log.sublog(), thresh=1e-15)
	except signal.ItimerError:
		out.log("Calculation Time Out!")
		return False
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
	return True





if __name__ == "__main__":
	resources, out  = parallel.resources(int(sys.argv[1])), output(log=textlog(echo=True))
	separation, n_ho = 5., 8
	for n_mol,max_states in [       [2,10],[3,10],[4,10],[5,10],[6,10],[7,10],[8,10],[9,10],[10,9],
	                         [11,9],[12,9],[13,8],[14,8],[15,7],[16,7],[17,7],[18,7],[19,6],[20,5],
	                         [21,5],[22,5],[23,5],[24,5],[25,3],[26,3],[27,3],[28,3],[29,3],[30,3]]:
		mol_obj = build(separation, n_mol, n_ho)
		analytical_value = do_analytical(mol_obj, out)
		if analytical_value is None:  out.log("+++++++++++++++++++++ %d ho %d mol           ++++++++++++++++++++++++++++++++++++++\nFAILED!!" %(n_ho, n_mol))
		else:
			n_state = 3
			while n_state <= max_states and main(mol_obj, [n_state for i in range(n_mol*n_ho)], out, analytical_value, resources):  n_state += 1
			if n_state < max_states:
				n_state += 1
				while n_state <= max_states:
					out.log("+++++++++++++++++++++ %d ho %d mol %d states  ++++++++++++++++++++++++++++++++++++++\nFAILED!!" %(n_ho, n_mol, n_state))
					n_state += 1
	out.log.write('primitive.out')
