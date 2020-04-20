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
import numpy as np
from qode.fermion_field.state import state, op_string, create, annihilate





def _digest_V_mat(V_mat_raw, num_spin_orb):
	V_mat = np.zeros( shape=(num_spin_orb, num_spin_orb, num_spin_orb, num_spin_orb) )
	for p in range(num_spin_orb):
		for q in range(num_spin_orb):
			for r in range(num_spin_orb):
				for s in range(num_spin_orb):				# Yikes!  V_mat_raw should come in as 4-index tensor!
					V_mat[p,q,r,s]  = V_mat_raw[p * num_spin_orb + q, s * num_spin_orb + r]
					V_mat[p,q,r,s] -= V_mat_raw[p * num_spin_orb + q, r * num_spin_orb + s]
					V_mat[p,q,r,s] /= 4
	return V_mat



# Some very thin wrappers that mimic functional programming style nested function definitions, because multiprocessing chokes on nested function defs.

class _one_elec_hamiltonian(object):
	#
	# H1 = \sum_pq h_pq p+ q
	# h_pq = Kinetic E integrals + Nuclear Attraction E integrals
	#
	def __init__(self,h_mat,num_spin_orb,configurations):
		self._arguments = (h_mat,num_spin_orb,configurations)
	def __call__(self,the_state):
		h_mat,num_spin_orb,configurations = self._arguments
		new_state = state(configurations)
		for p in range(num_spin_orb):
			for q in range(num_spin_orb):
				op = op_string( create(p), annihilate(q) )
				new_state.increment( op(the_state), h_mat[p,q] )
		return new_state

class _two_elec_hamiltonian(object):
	#
	# H2 = \sum_pqrs V_pqrs p+ q+ r s
	# V_pqrs = antisymmetrized 2 e- repulsion integrals
	#
	def __init__(self,V_mat_raw,num_spin_orb,configurations):
		V_mat = _digest_V_mat(V_mat_raw,num_spin_orb)
		self._arguments = (V_mat,num_spin_orb,configurations)
	def __call__(self,the_state):
		V_mat,num_spin_orb,configurations = self._arguments
		new_state = state(configurations)
		for p in range(num_spin_orb):
			for q in range(num_spin_orb):
				for r in range(num_spin_orb):
					for s in range(num_spin_orb):
						op = op_string( create(p), create(q), annihilate(r), annihilate(s) )
						new_state.increment( op(the_state), V_mat[p,q,r,s] )
		return new_state

class hamiltonian(object):
	def __init__(self,h_mat, V_mat_raw, configurations):
		num_spin_orb = h_mat.shape[0]
		H1 = _one_elec_hamiltonian(h_mat,    num_spin_orb,configurations)
		H2 = _two_elec_hamiltonian(V_mat_raw,num_spin_orb,configurations)
		self._arguments = (H1,H2,configurations)
	def __call__(self,the_state):
		H1,H2,configurations = self._arguments
		new_state = state(configurations)
		new_state.increment( H1(the_state) )
		new_state.increment( H2(the_state) )
		return new_state








"""\
def hamiltonian(state_obj, h_mat, V_mat_raw):
	occ_strings, _ = state_obj.value()
	H = hamiltonian_functional_programming_style(h_mat, V_mat_raw, occ_strings)
	return H(state_obj)


def h_wrapper(this_input_list):  # For CISD Code
	return hamiltonian(*this_input_list)


def H(state_obj, trivial_list_of_idx, trivial_t_amp_obj, trivial_occ_orbs, trivial_virt_orbs, trivial_offsets, h_mat, V_mat):
	return hamiltonian(state_obj, h_mat, V_mat)
"""









class inv_fock(object):
	def __init__(self, occ_orbs, vrt_orbs, fock_mat, configurations):
		fock_diag = [ fock_mat[i,i] for i in range(fock_mat.shape[0]) ]
		excitation_idx = self._get_excitation_idx(occ_orbs, vrt_orbs)
		orb_energy_diff_list = self._compute_orb_energy_diff(fock_diag, excitation_idx)
		self._arguments = (orb_energy_diff_list,configurations)
	def __call__(self,the_state):
		orb_energy_diff_list,configurations = self._arguments
		new_state = state(configurations)
		_,coefficients = the_state.value()
		coefficients = [c/d for c,d in zip(coefficients[1:],orb_energy_diff_list)]
		coefficients = [0.] + coefficients
		new_state.update_coeffs(coefficients)
		return new_state
	@staticmethod
	def _get_excitation_idx(occ_orbs, vrt_orbs):
		num_occ_orb = len(occ_orbs)
		num_vrt_orb = len(vrt_orbs)
		excitation = []
		for i in range(num_occ_orb):
			orb_i = occ_orbs[i]
			for a in range(num_vrt_orb):
				orb_a = vrt_orbs[a]
				excitation += [[orb_a,orb_i]]
		for i in range(num_occ_orb):
			orb_i = occ_orbs[i]
			for j in range(i+1, num_occ_orb):
				orb_j = occ_orbs[j]
				for a in range(num_vrt_orb):
					orb_a = vrt_orbs[a]
					for b in range(a+1, num_vrt_orb):
						orb_b = vrt_orbs[b]
						excitation += [[orb_a, orb_i, orb_b, orb_j]]
		return excitation
	@staticmethod
	def _compute_orb_energy_diff(fock_orb_energy, excitation_idx):
		energy_diff_list = []
		for item in excitation_idx:
			if len(item) == 2:
				energy_diff_list += [ fock_orb_energy[item[0]] - fock_orb_energy[item[1]] ]
			elif len(item) == 4:
				energy_diff_list += [ fock_orb_energy[item[0]] - fock_orb_energy[item[1]] + fock_orb_energy[item[2]] - fock_orb_energy[item[3]] ]
		return energy_diff_list


