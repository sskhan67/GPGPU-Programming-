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
from copy import deepcopy
import numpy as np
import qode
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



def _one_elec_hamiltonian(h_mat,num_spin_orb,configurations):
	#
	# H1 = \sum_pq h_pq p+ q
	# h_pq = Kinetic E integrals + Nuclear Attraction E integrals
	#
	def action(the_state):
		new_state = state(configurations)
		for p in range(num_spin_orb):
			for q in range(num_spin_orb):
				op = op_string( create(p), annihilate(q) )
				new_state.increment( op(the_state), h_mat[p,q] )
		return new_state
	return action

def _two_elec_hamiltonian(V_mat_raw,num_spin_orb,configurations):
	#
	# H2 = \sum_pqrs V_pqrs p+ q+ r s
	# V_pqrs = antisymmetrized 2 e- repulsion integrals
	#
	V_mat = _digest_V_mat(V_mat_raw,num_spin_orb)
	def action(the_state):
		new_state = state(configurations)
		for p in range(num_spin_orb):
			for q in range(num_spin_orb):
				for r in range(num_spin_orb):
					for s in range(num_spin_orb):
						op = op_string( create(p), create(q), annihilate(r), annihilate(s) )
						new_state.increment( op(the_state), V_mat[p,q,r,s] )
		return new_state
	return action

def hamiltonian(h_mat, V_mat_raw, configurations):
	num_spin_orb = h_mat.shape[0]
	H1 = _one_elec_hamiltonian(h_mat,    num_spin_orb,configurations)
	H2 = _two_elec_hamiltonian(V_mat_raw,num_spin_orb,configurations)
	def action(the_state):
		new_state = state(configurations)
		new_state.increment( H1(the_state) )
		new_state.increment( H2(the_state) )
		return new_state
	return action





if __name__ == "__main__":
	occ_orbs = [0, 4]
	vrt_orbs = [1,2,3, 5,6,7]
	configurations = qode.fermion_field.occ_strings.CISD(occ_orbs,vrt_orbs)
	#
	h_mat = np.load('data/h_mat.npy')
	V_mat = np.load('data/V_mat.npy')
	H1 = _one_elec_hamiltonian(h_mat,8,configurations)
	H2 = _two_elec_hamiltonian(V_mat,8,configurations)
	H = hamiltonian(h_mat,V_mat,configurations)
	#
	phi0 = state(configurations,0)
	#
	print("|Phi0> = ")
	phi0.print()
	psi1 = H1(phi0)
	print("H1|Phi0> = ")
	psi1.print()
	"""\
	"""
	psi2 = H2(phi0)
	print("H2|Phi0> = ")
	psi2.print()
	psi = deepcopy(psi1)
	psi.increment(psi2)
	print("H1|Phi0> + H2|Phi0> = ")
	psi.print()
	print("H|Phi0> = ")
	H(phi0).print()
