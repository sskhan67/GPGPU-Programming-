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
from qode.fermion_field.state import state, op_string, create, annihilate, dot
from qode.fermion_field.reorder import *


class operator(object):
	def __init__(self,configurations,ic_range,id_range,ac_range,ad_range,C0=0.):
		self.C0 = C0
		self.Ex = None
		self.Dx = None
		self.Fo = None
		self.Fv = None
		self.ic_range = ic_range
		self.ac_range = ac_range
		self.id_range = id_range
		self.ad_range = ad_range
		self.configurations = configurations
	def __call__(self,the_state):
		new_state = state(self.configurations)
		new_state.increment(the_state,self.C0)
		if self.Ex is not None:
			for a in self.ac_range:
				for i in self.id_range:
					op = op_string( create(a), annihilate(i) )
					new_state.increment( op(self.Ex[a][i](the_state)) )
		if self.Dx is not None:
			for a in self.ad_range:
				for i in self.ic_range:
					op = op_string( annihilate(a), create(i) )
					new_state.increment( op(self.Dx[a][i](the_state)) )
		if self.Fo is not None:
			for i in self.id_range:
				for j in self.ic_range:
					op = op_string( annihilate(i), create(j) )
					new_state.increment( op(self.Fo[i][j](the_state)) )
		if self.Fv is not None:
			for a in self.ac_range:
				for b in self.ad_range:
					op = op_string( create(a), annihilate(b) )
					new_state.increment( op(self.Fv[a][b](the_state)) )
		return new_state
	def add_this_to(self,other):
		other.C0 += self.C0
		if self.Ex is not None:
			target = other._get_Ex_block()
			for a in self.ac_range:
				for i in self.id_range:
					self.Ex[a][i].add_this_to(target[a][i])
		if self.Dx is not None:
			target = other._get_Dx_block()
			for a in self.ad_range:
				for i in self.ic_range:
					self.Dx[a][i].add_this_to(target[a][i])
		if self.Fo is not None:
			target = other._get_Fo_block()
			for i in self.id_range:
				for j in self.ic_range:
					self.Fo[i][j].add_this_to(target[i][j])
		if self.Fv is not None:
			target = other._get_Fv_block()
			for a in self.ac_range:
				for b in self.ad_range:
					self.Fv[a][b].add_this_to(target[a][b])
	def increment(self,other):  other.add_this_to(self)
	def _get_Ex_block(self):
		if self.Ex is not None:  return self.Ex
		else:
			self.Ex = {}
			for na,a in enumerate(self.ac_range):
				row = {}
				for ni,i in enumerate(self.id_range):  row[i] = operator(self.configurations,self.ic_range,self.id_range[ni+1:],self.ac_range[na+1:],self.ad_range)
				self.Ex[a] = row
			return self.Ex
	def _get_Dx_block(self):
		if self.Dx is not None:  return self.Dx
		else:
			self.Dx = {}
			for na,a in enumerate(self.ad_range):
				row = {}
				for ni,i in enumerate(self.ic_range):  row[i] = operator(self.configurations,self.ic_range[ni+1:],None,None,self.ad_range[na+1:])
				self.Dx[a] = row
			return self.Dx
	def _get_Fo_block(self):
		if self.Fo is not None:  return self.Fo
		else:
			self.Fo = {}
			for ni,i in enumerate(self.id_range):
				row = {}
				for nj,j in enumerate(self.ic_range):  row[j] = operator(self.configurations,self.ic_range[nj+1:],self.id_range[ni+1:],None,self.ad_range)
				self.Fo[i] = row
			return self.Fo
	def _get_Fv_block(self):
		if self.Fv is not None:  return self.Fv
		else:
			self.Fv = {}
			for na,a in enumerate(self.ac_range):
				row = {}
				for nb,b in enumerate(self.ad_range):  row[b] = operator(self.configurations,self.ic_range,None,self.ac_range[na+1:],self.ad_range[nb+1:])
				self.Fv[a] = row
			return self.Fv
	def add_term(self,term):
		try:
			self.C0 += term[0]
		except:
			if   term[0].is_excitation:
				a,i = term[0].E_transition
				if a not in self.ac_range or i not in self.id_range:  raise an_error
				Ex = self._get_Ex_block()
				Ex[a][i].add_term(term[1:])
			elif term[0].is_deexcitation:
				a,i = term[0].D_transition
				if a not in self.ad_range or i not in self.ic_range:  raise an_error
				Dx = self._get_Dx_block()
				Dx[a][i].add_term(term[1:])
			elif term[0].is_occ_rearrange:
				i,j = term[0].Fo_transition
				if i not in self.id_range or j not in self.ic_range:  raise an_error
				Fo = self._get_Fo_block()
				Fo[i][j].add_term(term[1:])
			elif term[0].is_vrt_rearrange:
				a,b = term[0].Fv_transition
				if a not in self.ac_range or b not in self.ad_range:  raise an_error
				Fv = self._get_Fv_block()
				Fv[a][b].add_term(term[1:])



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



def hamiltonian(h_mat, V_mat_raw, configurations, occ_orbs, vrt_orbs):
	num_orb = len(occ_orbs) + len(vrt_orbs)
	V_mat = _digest_V_mat(V_mat_raw,num_orb)
	H = operator(configurations,occ_orbs,occ_orbs,vrt_orbs,vrt_orbs)
	reorder_H1 = reorder_CpAq(vrt_orbs)
	reorder_H2 = reorder_CpCqArAs(vrt_orbs)
	for p in range(num_orb):
		for q in range(num_orb):
			terms = reorder_H1(p,q,h_mat[p,q])
			for term in terms:
				phase = term.pop(0)	# this should be handled . . .
				term.append(phase)	# . . . at a lower level
				H.add_term(term)
	for p in range(num_orb):
		for q in range(num_orb):
			for r in range(num_orb):
				for s in range(num_orb):
					terms = reorder_H2(p,q,r,s,V_mat[p,q,r,s])
					for term in terms:
						phase = term.pop(0)	# this should be handled . . .
						term.append(phase)	# . . . at a lower level
						H.add_term(term)
	return H



def SD_excitations(reference, the_state, configurations, occ_orbs, vrt_orbs):
	X = operator(configurations,None,occ_orbs,vrt_orbs,None)
	for a in vrt_orbs:
		for i in occ_orbs:
			amplitude = dot(op_string(create(a),annihilate(i))(reference),the_state)
			X.add_term([excitation(a,i),amplitude])
	for idxa,a in enumerate(vrt_orbs):
		for idxi,i in enumerate(occ_orbs):
			for b in vrt_orbs[idxa+1:]:
				for j in occ_orbs[idxi+1:]:
					amplitude = dot(op_string(create(a),annihilate(i),create(b),annihilate(j))(reference),the_state)
					X.add_term([excitation(a,i),excitation(b,j),amplitude])
	return X
