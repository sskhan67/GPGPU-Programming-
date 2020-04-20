#    (C) Copyright 2018, 2019 Anthony D. Dutoi and Yuhong Liu
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



def dist_mat(separation, N_frag):
	txt_matrix = ""		# The matrix is text based because the distances are primarily for labeling files (though they can also be parsed to floats)
	prepend = ""
	for n in range(N_frag):
		row = ""
		for m in range(n+1):
			row += "    . "
		sep = separation
		for m in range(n+1,N_frag):
			row += " " + "{:.1f}".format(sep).rjust(5)
			sep += separation
		txt_matrix += prepend+row
		prepend = "\n"
	txt_matrix = [line.split() for line in txt_matrix.split("\n")]	# sort of silly, but the above prints nicely, in case I want to resurrect that
	return txt_matrix

def E_nuc(dist_mat, Z):
	E = 0
	for dist_row in dist_mat:
		for dist in dist_row:
			if dist!=".":
				R = float(dist)
				E += Z**2 * (0.52917721092 / R)
	return E



def _ket_loop(Hmat, bra_states, I, ket_states_now, J, N_frag, states_per_frag, H):
	n = len(ket_states_now) + 1
	ket_states_next = list(ket_states_now) + [None]
	if n==N_frag:
		monomer_Hamiltonians, dimer_Couplings = H
		ket_states = ket_states_next
		for i in range(states_per_frag):
			ket_states[-1] = i
			transitions = []
			for m,(bra,ket) in enumerate(zip(bra_states,ket_states)):
				if bra!=ket:  transitions += [(m,bra,ket)]
			if len(transitions)==0:
				for M in range(0,N_frag):
					Hmat[I,J] += monomer_Hamiltonians[M][bra_states[M], ket_states[M]]
					for N in range(M+1,N_frag):
						Hmat[I,J] += dimer_Couplings[M][N][bra_states[M]*states_per_frag+bra_states[N], ket_states[M]*states_per_frag+ket_states[N]]
			if len(transitions)==1:
				(M, Mb, Mk), = transitions
				Hmat[I,J] += monomer_Hamiltonians[M][Mb,Mk]
				for N in range(0,M):
					Hmat[I,J] += dimer_Couplings[N][M][bra_states[N]*states_per_frag+Mb, ket_states[N]*states_per_frag+Mk]
				for N in range(M+1,N_frag):
					Hmat[I,J] += dimer_Couplings[M][N][Mb*states_per_frag+bra_states[N], Mk*states_per_frag+ket_states[N]]
			if len(transitions)==2:
				(M, Mb, Mk),(N, Nb, Nk) = transitions
				Hmat[I,J] += dimer_Couplings[M][N][Mb*states_per_frag+Nb, Mk*states_per_frag+Nk]
			J += 1
	else:
		stride = states_per_frag**(N_frag-n)
		for i in range(states_per_frag):
			ket_states_next[-1] = i
			_ket_loop(Hmat, bra_states, I, ket_states_next, J, N_frag, states_per_frag, H)
			J += stride

def _bra_loop(Hmat, bra_states_now, I, N_frag, states_per_frag, H):
	n = len(bra_states_now) + 1
	bra_states_next = list(bra_states_now) + [None]
	if n==N_frag:
		bra_states = bra_states_next
		for i in range(states_per_frag):
			bra_states[-1] = i
			_ket_loop(Hmat, bra_states, I, [], 0, N_frag, states_per_frag, H)
			I += 1
	else:
		stride = states_per_frag**(N_frag-n)
		for i in range(states_per_frag):
			bra_states_next[-1] = i
			_bra_loop(Hmat, bra_states_next, I, N_frag, states_per_frag, H)
			I += stride

# Built on the assumption that states_per_frag is the same for all fragments
def braket_loops(Hmat, N_frag, states_per_frag, H):
	_bra_loop(Hmat, [], 0, N_frag, states_per_frag, H)
