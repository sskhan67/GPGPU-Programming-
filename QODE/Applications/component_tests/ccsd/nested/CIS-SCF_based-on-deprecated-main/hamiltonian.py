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
from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
from qode.many_body.nested_operator import operator, mask, multiply



def digest_V_mat(V_mat_raw, num_spin_orb):
	V_mat = np.zeros( shape=(num_spin_orb, num_spin_orb, num_spin_orb, num_spin_orb) )
	for p in range(num_spin_orb):
		for q in range(num_spin_orb):
			for r in range(num_spin_orb):
				for s in range(num_spin_orb):				# Yikes!  V_mat_raw should come in as 4-index tensor!
					V_mat[p,q,r,s]  = V_mat_raw[p * num_spin_orb + q, s * num_spin_orb + r]
					V_mat[p,q,r,s] -= V_mat_raw[p * num_spin_orb + q, r * num_spin_orb + s]
					V_mat[p,q,r,s] /= 4
	return V_mat



def hamiltonian(E_nuc, h_mat, V_mat, fluctuations):
	H = operator(fluctuations)
	occ_orbs, vrt_orbs = fluctuations.occ_orbs, fluctuations.vrt_orbs
	num_orb = len(occ_orbs) + len(vrt_orbs)
	reorder_H1 = reorder_CpAq(vrt_orbs)
	reorder_H2 = reorder_CpCqArAs(vrt_orbs)
	H.add_term([E_nuc])
	for p in range(num_orb):
		for q in range(num_orb):
			terms = reorder_H1(p,q,h_mat[p,q])
			for term in terms:  H.add_term(term)
	for p in range(num_orb):
		for q in range(num_orb):
			for r in range(num_orb):
				for s in range(num_orb):
					terms = reorder_H2(p,q,r,s,V_mat[p,q,r,s])
					for term in terms:  H.add_term(term)
	return H



class PT_step_conditioner(object):
	def __init__(self, energy_diffs, fluctuations):
		self.fluctuations = fluctuations
		self.energy_diffs = energy_diffs
	def __call__(self, Omega):
		new = operator(self.fluctuations)
		self._apply_denominators(new, Omega)
		return new
	def _apply_denominators(self, new, old, shift=0):		# Very dirty, looking at another object's privates
		if old._Ex is not None:
			new._get_Ex_block()
			for index in self.fluctuations.all_transition_indices().Excitations:
				denom = self.energy_diffs[index] + shift
				new._Ex[index]._C[0] = old._Ex[index]._C[0] / denom			# Dividing ... 
				self._apply_denominators(new._Ex[index], old._Ex[index], shift=denom)	# ... and then recurring misses reference, which is desired.



excite_mask = mask( [False,[True,False,False,False],False,False] )

class const_excite_projection_class(object):
	""" wraps a mask so that it behaves like a callable projection operation on an operator, returning a new operator, for use with operator spaces """
	def __init__(self):
		self.mask = mask( [True,[True,False,False,False],False,False] )
	def __call__(self, op):
		print("... and projected")
		return self.mask(op).copy()

const_excite_projection = const_excite_projection_class()

class super_operator_wrapper(object):
	""" wraps an operator so that it behaves like a callable operation on another operator, returning a new operator, for use with operator spaces """
	def __init__(self, op, resources):
		self.op        = op
		self.resources = resources
	def __call__(self, other):
		print("H acted ...")
		return multiply(self.op, other, resources=self.resources)
