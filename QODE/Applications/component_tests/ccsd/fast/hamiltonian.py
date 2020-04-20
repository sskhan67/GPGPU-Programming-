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
from qode.util import parallel
from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
from qode.many_body.nested_operator import operator as nested_operator
from qode.many_body.nested_operator import aggregator
from qode.many_body.fast_operator import operator
from qode.many_body.fast_operator.generalutils import load_pickle_into_blocks


from qode.fermion_field.reorder import excitation    as Ex
from qode.fermion_field.reorder import deexcitation  as Dx
from qode.fermion_field.reorder import occ_rearrange as Fo
from qode.fermion_field.reorder import vrt_rearrange as Fv

np.set_printoptions(threshold=np.nan,precision=6)


def op_get_idx(trans_obj, No, Nv):
	if isinstance(trans_obj, Ex):
		a, i = trans_obj.E_transition
		return [i,a], [No,Nv], 'Ex'
	elif isinstance(trans_obj, Dx):
		a, i = trans_obj.D_transition
		return [a,i], [Nv,No], 'Dx'
	elif isinstance(trans_obj, Fo):
		i, j = trans_obj.Fo_transition
		return [i,j], [No,No], 'Fo'
	elif isinstance(trans_obj, Fv):
		a, b = trans_obj.Fv_transition
		return [b,a], [Nv,Nv], 'Fv'
	else:
		raise AssertionError


def parse_term(term, No, Nv):
	# Take an operator term, clean it up to my storage format 
	# (subscript looped first, superscript looped later)
	if len(term) == 1:
		if term[0] == 0:
			return None, None, None
	elif len(term) > 1:
		sign = term[-1]
		all_idx = []
		dims    = []
		block   = ''
		for item in term[:-1]:
			idx, dim, trans = op_get_idx(item, No, Nv)
			all_idx += idx
			dims    += dim
			block   += trans
		return block, all_idx, sign
	else:
		raise AssertionError


def myprint(mat):
	dim = 6
	output_buf = ""
	for i in range(dim**2):
		for j in range(dim**2):
			if mat[i,j] != 0.0:
				output_buf += "%.6e\t" %(mat[i,j])
			else:
				output_buf += "0.\t"
		output_buf += '\n'
	print(output_buf)

def myprint_dig(mat):
	dim = 6
	output_buf = ""
	for i in range(dim):
		for j in range(dim):
			for k in range(dim):
				for l in range(dim):
					if mat[i,j,k,l] != 0.0:
						output_buf += "%.6e\t" %(mat[i,j,k,l])
					else:
						output_buf += "0.\t"
			output_buf += '\n'
	print(output_buf)


def _digest_V_mat(V_mat_raw, num_spin_orb):
	# print("RAW V...")
	# myprint(V_mat_raw)
	# print("-----------------------------------------------------------------------------------")
	V_mat = np.zeros( shape=(num_spin_orb, num_spin_orb, num_spin_orb, num_spin_orb) )
	for p in range(num_spin_orb):
		for q in range(num_spin_orb):
			for r in range(num_spin_orb):
				for s in range(num_spin_orb):				# Yikes!  V_mat_raw should come in as 4-index tensor!
					V_mat[p,q,r,s]  = V_mat_raw[p * num_spin_orb + q, s * num_spin_orb + r]
					V_mat[p,q,r,s] -= V_mat_raw[p * num_spin_orb + q, r * num_spin_orb + s]
					V_mat[p,q,r,s] /= 4
	# print("NEW V...")
	# myprint_dig(V_mat)
	# print("-----------------------------------------------------------------------------------")
	return V_mat



def partial_h(p_range, h_mat, fluctuations):
	#print("starting h",p_range)
	partial = nested_operator(fluctuations)
	occ_orbs, vrt_orbs = fluctuations.occ_orbs, fluctuations.vrt_orbs
	num_orb = len(occ_orbs) + len(vrt_orbs)
	reorder_H1 = reorder_CpAq(vrt_orbs)
	#print("looping h",p_range)
	for p in p_range:
		for q in range(num_orb):
			terms = reorder_H1(p,q,h_mat[p,q])
			#
			# for term in terms:
			# 	print('V[',p,q, '] -->', parse_term( term, len(occ_orbs), len(vrt_orbs) ) )
			for term in terms:  partial.add_term(term)
	#print("done h", p_range)
	return partial

def partial_V(p_range, V_mat, fluctuations):
	#print("starting V",p_range)
	partial = nested_operator(fluctuations)
	occ_orbs, vrt_orbs = fluctuations.occ_orbs, fluctuations.vrt_orbs
	num_orb = len(occ_orbs) + len(vrt_orbs)
	reorder_H2 = reorder_CpCqArAs(vrt_orbs)
	#print("looping V",p_range)
	for p in p_range:
		for q in range(num_orb):
			for r in range(num_orb):
				for s in range(num_orb):
					terms = reorder_H2(p,q,r,s,V_mat[p,q,r,s])
					#
					# for term in terms:
					# 	print( 'V[', p,q,r,s, '] -->', parse_term( term, len(occ_orbs), len(vrt_orbs) ) )
					for term in terms:  partial.add_term(term)
	#print("done V", p_range)
	return partial

def hamiltonian(E_nuc, h_mat, V_mat_raw, fluctuations, resources=parallel.resources(1)):
	occ_orbs, vrt_orbs = fluctuations.occ_orbs, fluctuations.vrt_orbs
	num_orb = len(occ_orbs) + len(vrt_orbs)
	V_mat = _digest_V_mat(V_mat_raw,num_orb)
	H = nested_operator(fluctuations)
	H.add_term([E_nuc])
	n_cores = resources.n_cores
	#if n_cores==1:
	if True:			# temporary hack to make debugging more sane
		for p_range in [range(num_orb)]:
			H.increment(partial_h(p_range, h_mat, fluctuations))
		for p_range in [range(num_orb)]:
			H.increment(partial_V(p_range, V_mat, fluctuations))
	else:
		if n_cores<=num_orb:
			chunk = num_orb // n_cores	# chunk will always be a little too small
			p_ranges = []
			begin = 0
			for core in range(num_cores-1):
				end = begin + chunk
				p_ranges += [range(begin,end)]
				begin = end
			p_ranges += [range(begin,num_orb)]
		else:
			n_cores = num_orb
			p_ranges = [ range(p,p+1) for p in range(num_orb) ]
		p_ranges = [(p_range,) for p_range in p_ranges]
		H = parallel.aggregate(partial_h, (h_mat,fluctuations), aggregator(H), n_cores).run(H, p_ranges)
		H = parallel.aggregate(partial_V, (V_mat,fluctuations), aggregator(H), n_cores).run(H, p_ranges)
	fast_H = operator()
	load_pickle_into_blocks(H.dump(), fast_H, occ_orbs, vrt_orbs)
	return fast_H
