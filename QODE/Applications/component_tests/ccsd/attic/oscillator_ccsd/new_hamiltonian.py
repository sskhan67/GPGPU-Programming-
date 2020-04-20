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
import sys, os, time
import numpy as np
from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
from qode.many_body.nested_operator import operator as nested_operator
from qode.many_body.fast_operator.generalutils  import print_array, add_value_to_block, load_pickle_into_blocks
from qode.many_body.fast_operator.space_traits  import add_cluster_ops
from qode.many_body.fast_operator.fast_operator import cc_operator
from qode                                       import fermion_field
from qode.many_body.fast_operator       import operator, space_traits, BakerCampbellHausdorff
import copy
from copy import deepcopy
import pickle
from qode.many_body                     import CCSD
from qode.util import parallel, output, textlog, indent

def init_H_blocks_zero(cc_op_obj, No, Nv):
	cc_op_obj.K      = np.zeros(1)        # This is zero
	cc_op_obj.Ex     = np.zeros((No*Nv))  # This is zero
	cc_op_obj.Fo     = np.zeros((No*No))
	cc_op_obj.Fv     = np.zeros((Nv*Nv))
	cc_op_obj.Dx     = np.zeros((Nv*No))  # This is zero
	cc_op_obj.ExEx   = np.zeros(((No*Nv)**2))
	cc_op_obj.ExFo   = np.zeros((Nv*No**3))
	cc_op_obj.ExFv   = np.zeros((No*Nv**3))
	cc_op_obj.ExDx   = np.zeros(((No*Nv)**2))
	cc_op_obj.FoFo   = np.zeros((No**4))
	cc_op_obj.FvFv   = np.zeros((Nv**4))
	cc_op_obj.FoDx   = np.zeros((Nv*No**3))
	cc_op_obj.FvDx   = np.zeros((No*Nv**3))
	cc_op_obj.DxDx   = np.zeros(((No*Nv)**2))	
	return cc_op_obj



def _digest_V_mat(V_mat_raw, num_spin_orb, rec_num_states):
	from_mol = []
	for i in range(len(rec_num_states)):
		from_mol += [i for n in range(rec_num_states[i])]
	print("FROM MOL =",from_mol)
	mol_loc = []
	for n in rec_num_states:
		mol_loc += [i for i in range(n)]
	print("MOL LOC=",mol_loc)

	V_mat_full = np.zeros((num_spin_orb*num_spin_orb,num_spin_orb*num_spin_orb))
	for p in range(num_spin_orb):
		for q in range(num_spin_orb):
			for r in range(num_spin_orb):
				for s in range(num_spin_orb):
					if from_mol[p] == from_mol[r] and from_mol[q] == from_mol[s]:
						if from_mol[p] < from_mol[q]:
							V_mat_full[ p*num_spin_orb+q , r*num_spin_orb+s ] = V_mat_raw(from_mol[p], from_mol[q], mol_loc[p], mol_loc[q], mol_loc[s], mol_loc[r])
						elif from_mol[p] > from_mol[q]:
							V_mat_full[ p*num_spin_orb+q , r*num_spin_orb+s ] = V_mat_raw(from_mol[q], from_mol[p], mol_loc[s], mol_loc[r], mol_loc[p], mol_loc[q])




	V_mat = np.zeros( shape=(num_spin_orb, num_spin_orb, num_spin_orb, num_spin_orb) )
	for p in range(num_spin_orb):
		for q in range(num_spin_orb):
			for r in range(num_spin_orb):
				for s in range(num_spin_orb):				# Yikes!  V_mat_raw should come in as 4-index tensor!
					V_mat[p,q,r,s]  = V_mat_full[p * num_spin_orb + q, s * num_spin_orb + r]
					V_mat[p,q,r,s] -= V_mat_full[p * num_spin_orb + q, r * num_spin_orb + s]
					V_mat[p,q,r,s] /= 4						
	return V_mat

def hamiltonian(h_mat, V_mat_raw, fluctuations, rec_num_states):
	H = nested_operator(fluctuations)
	occ_orbs, vrt_orbs = fluctuations.occ_orbs, fluctuations.vrt_orbs
	num_orb = len(occ_orbs) + len(vrt_orbs)
	V_mat = _digest_V_mat(V_mat_raw,num_orb,rec_num_states)
	reorder_H1 = reorder_CpAq(vrt_orbs)
	reorder_H2 = reorder_CpCqArAs(vrt_orbs)
	# H.add_term([E_nuc])
	for p in range(num_orb):
		for q in range(num_orb):
			terms = reorder_H1(p,q,h_mat(p,q))  # my h_mat has a __call__ function takes (p,q)
			for term in terms:  H.add_term(term)
	for p in range(num_orb):
		for q in range(num_orb):
			for r in range(num_orb):
				for s in range(num_orb):
					terms = reorder_H2(p,q,r,s,V_mat[p,q,r,s])
					for term in terms:  H.add_term(term)
	fast_H = operator()
	load_pickle_into_blocks(H.dump(), fast_H, occ_orbs, vrt_orbs )
	print_array(fast_H)
	return fast_H

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

class orbital_energy(object):
	"""Class for modifying Omega operators with orbital energy differences"""
	def __init__(self, all_mol_energy, occ_orbs, vrt_orbs):
		self.orbital_energy = all_mol_energy
		self.No = len(occ_orbs)
		self.Nv = len(vrt_orbs)
		self.occ_orbs = occ_orbs
		self.vrt_orbs = vrt_orbs
	def __call__(self, omega_obj):
		inverted_omega = deepcopy(omega_obj)
		print_array(inverted_omega)
		inverted_omega.K[0] = 0.0
		for i in range(self.No):
			for a in range(self.Nv):
				inverted_omega.Ex[i*self.Nv + a] /= (self.orbital_energy[self.vrt_orbs[a]] -  self.orbital_energy[self.occ_orbs[i]]) #* -1.0
		for i in range(self.No):
			for a in range(self.Nv):
				for j in range(i):
					for b in range(a):
						inverted_omega.ExEx[i*self.Nv*self.No*self.Nv + a*self.No*self.Nv + j* self.Nv + b] /= \
						        (self.orbital_energy[self.vrt_orbs[a]] + self.orbital_energy[self.vrt_orbs[b]] -  self.orbital_energy[self.occ_orbs[i]] - self.orbital_energy[self.occ_orbs[j]]) #*-1.0
		return inverted_omega


if __name__ == "__main__":
	from Applications.component_tests.oscillator_ccsd.old_qode.generate_ho_config import diag_rec_h0_mat, diag_rec_v_matrix, generate_main, get_lowest_analytical, get_rec_obj
	# Oscillator System Imports
	from qode.coupled_oscillators            import oscillator_system as osys
	from qode.coupled_oscillators            import hamiltonian_class
	from qode.coupled_oscillators.analytical import analytical_solution



	rec_num_states = [3,3,3,3]
	
	ho1 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho2 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	

	ho3 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho4 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )
	
	ho5 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho6 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )

	ho7 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 1.0 ) )
	ho8 = osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = 2.0 ) )

	
	coupling_mat = [[0.0, 0.3],\
	                [0.3, 0.0]]

	# coupling_mat = 
	
	mol1 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho1,ho2] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	
	mol2 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho3,ho4] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	
	mol3 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho5,ho6] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )
	
	mol4 = osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( [ho7,ho8] ),
                                                            coupling       = copy.deepcopy( coupling_mat )  ) )


	# ext_k_mat = [[0.0,0.1],[0.1,0.0]]
	ext_k_mat =    [[0.0, 0.1, 0.1, 0.1],\
                    [0.0, 0.0, 0.1, 0.1],\
                    [0.0, 0.0, 0.0, 0.1],\
                    [0.0, 0.0, 0.0, 0.0]]
	
	mol_obj = osys.OscillatorSystem( osys.group_of_oscillators( recursive_list =  [mol1, mol2, mol3, mol4],
                                                                coupling       =  ext_k_mat ) )
	
	rec_obj = get_rec_obj( mol_obj, rec_num_states )
	
	
	list_mol_eigen_values = rec_obj.get_sorted_eigen_energy()
	v_mat_list = rec_obj.get_pair_dipole_mat()
	

	# Building h_mat and V_mat objects
	h_mat = diag_rec_h0_mat(list_mol_eigen_values)
	V_mat = diag_rec_v_matrix(rec_num_states, v_mat_list)





	No = len(rec_num_states)
	Nv = 0
	for item in rec_num_states:
		Nv += (item - 1)
	#
	occ_orbs = [0]
	vrt_orbs = [ occ_orbs[-1]+1+i for i in range(rec_num_states[0]-1) ]
	for n in rec_num_states[1:]:
		occ_orbs += [ vrt_orbs[-1]+1 ]
		vrt_orbs += [ occ_orbs[-1]+1+i for i in range(n-1) ]
	print("OCC =", No, occ_orbs)
	print("VRT =", Nv, vrt_orbs)

	print("H original dimension =", h_mat.get_dimension())
	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)




	print("ANALYTICAL SOLUTION =", get_analytical_solution(mol_obj))
	H = hamiltonian(h_mat, V_mat, fluctuations, rec_num_states)
	# pickle.dump(H, open("3mol.pickled",'wb'))

	#
	# 3 molecules E = 3.3816016755152276
	#
	out = output(log=textlog(echo=True))
	out = output(log=out.log.sublog())
	out.log("Parsing input ...")

	out.log(indent("occupied orbital indices = ", occ_orbs))
	out.log(indent("virtual  orbital indices = ", vrt_orbs))

	out.log("Partioning fluctuations (excitations, flat, deexcitations) ...")
	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)

	out.log("Building normal-ordered H ...")
	# H = pickle.load( open("3mol.pickled", "rb") )

	out.log("Getting orbital energy preconditioner ...")
	all_mol_energy = []
	for item in list_mol_eigen_values:
		all_mol_energy += item
	# print("ALL ENERGY =", all_mol_energy)
	invF = orbital_energy(all_mol_energy, occ_orbs, vrt_orbs)	# DANGER hardcoded number of e-

	out.log("Initializing T ...")
	T = operator()			# Initialized to zero

	resources = []
	BCH = BakerCampbellHausdorff(len(occ_orbs), len(vrt_orbs), resources)
	out.log("Baker-Campbell-Hausdorff module initialized.")

	out.log("Perform coupled-cluster iterations.")
	t_ccsd_start = time.time()
	E = CCSD(H, T, invF, BCH, space_traits, out.log.sublog())
	t_ccsd_end = time.time()

	out.log("CCSD Energy =", E, 'Hartree')
	out.log("CCSD TIME {}".format(t_ccsd_end - t_ccsd_start))