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
from qode.util import parallel
from qode.fermion_field.state import *



def increment(target,delta):
	target.increment(delta)
	return target

def subcall(op,coeff,the_state):
	delta = op | the_state
	delta.scale(coeff)
	return delta


class H(object):
	def __init__(self, h, V, occ_strings):
		self.occ_strings = occ_strings
		n_spatial_orb = h.shape[0]
		spatial_indices = range(n_spatial_orb)
		self.terms = []
		for p in spatial_indices:
			for p_spin in ("dn","up"):
				for q in spatial_indices:
					for q_spin in ("dn","up"):
						if p_spin==q_spin:
							offset = 0 if p_spin=="dn" else 1
							pq = op_string( create(2*p+offset), annihilate(2*q+offset) )
							self.terms += [( pq, h[p,q] )]
		for p in spatial_indices:
			for p_spin in ("dn","up"):
				for q in spatial_indices:
					for q_spin in ("dn","up"):
						for r in spatial_indices:
							for r_spin in ("dn","up"):
								for s in spatial_indices:
									for s_spin in ("dn","up"):
										if p_spin==s_spin and q_spin==r_spin:
											offset1 = 0 if p_spin=="dn" else 1
											offset2 = 0 if q_spin=="dn" else 1
											pqrs = op_string( create(2*p+offset1), create(2*q+offset2), annihilate(2*r+offset2), annihilate(2*s+offset1) )
											self.terms += [( pqrs, V[p,s,q,r]/2 )]
										#if p_spin==r_spin and q_spin==s_spin:
										#	offset1 = 0 if p_spin=="dn" else 1
										#	offset2 = 0 if q_spin=="dn" else 1
										#	pqrs = op_string( create(2*p+offset1), create(2*q+offset2), annihilate(2*r+offset1), annihilate(2*s+offset2) )
										#	self.terms += [( pqrs, -V[p,r,q,s]/4 )]
	def __call__(self, the_state):
		result = state(self.occ_strings)
		#for op,coeff in self.terms:  result.increment( op(the_state), coeff )
		n_cores = 1
		if n_cores==1:
			for op,coeff in self.terms:
				result = increment(result, subcall(op,coeff,the_state))
		else:
			result = parallel.aggregate(subcall, (the_state,), increment, n_cores).run(result, self.terms)
		return result
