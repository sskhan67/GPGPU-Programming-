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
from qode.util.parallel import parallelize_task
from state import *



def couple_transition_to_other_molecule(M, m_f, m_i, N, num_states_N, transition_rhos, basis_ranges, V, nuc_mats):
	rho_M = transition_rhos[M][m_f][m_i]

	terms = []
	for n_f in range(num_states_N):
		for n_i in range(num_states_N):
			rho_N = transition_rhos[N][n_f][n_i]

			V_ = 0
			for i_p,p in enumerate(range(*basis_ranges[M])):
				for i_s,s in enumerate(range(*basis_ranges[M])):
					for i_q,q in enumerate(range(*basis_ranges[N])):
						for i_r,r in enumerate(range(*basis_ranges[N])):
							V_ += rho_M[i_p,i_s] * V[p,s,q,r] * rho_N[i_q,i_r]
			if m_f==m_i:
				for i_q,q in enumerate(range(*basis_ranges[N])):
					for i_r,r in enumerate(range(*basis_ranges[N])):
						V_ += nuc_mats[N][i_q,i_r] * rho_N[i_q,i_r]
			if n_f==n_i:
				for i_p,p in enumerate(range(*basis_ranges[M])):
					for i_s,s in enumerate(range(*basis_ranges[M])):
						V_ += rho_M[i_p,i_s] * nuc_mats[M][i_p,i_s]
				
			op = op_string( transition(M,m_i,m_f) , transition(N,n_i,n_f) )
			terms += [( op, V_ )]

	return terms



class H(object):
	def __init__(self, eigen_energies, transition_rhos, basis_ranges, V, nuc_mats, excite_strings):
		self.excite_strings = excite_strings
		num_molecules = len(eigen_energies)

		self.terms = []

		for M in range(num_molecules):
			num_states = len(eigen_energies[M])

			for m in range(num_states):
				h_ = eigen_energies[M][m]
				op = op_string( transition(M,m,m) )
				self.terms += [( op, h_ )]

		inputs = []
		for M in range(num_molecules):
			num_states_M = len(eigen_energies[M])
			for N in range(0,M):
				num_states_N = len(eigen_energies[N])

				for m_f in range(num_states_M):
					for m_i in range(num_states_M):
						inputs += [(M, m_f, m_i, N, num_states_N, transition_rhos, basis_ranges, V, nuc_mats)]

		#for arguments in inputs:
		#	self.terms += couple_transition_to_other_molecule(*arguments)
		nested_terms = parallelize_task(couple_transition_to_other_molecule, inputs, n_workers=25)
		for batch in nested_terms:  self.terms += batch

		print("Hamiltonian loaded")

	def __call__(self, the_state):
		result = state(self.excite_strings)
		for op,coeff in self.terms:  result.increment( op|the_state, coeff )
		return result
