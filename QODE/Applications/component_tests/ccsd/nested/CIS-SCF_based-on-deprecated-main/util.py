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
import numpy
from qode.math import linear_inner_product_space, numpy_space, sqrt
from qode.fermion_field import OV_partitioning



def symmetric_ON(vecs):
	n = len(vecs)
	S = numpy.zeros((n,n))

	for Np,Vp in enumerate(vecs):
		for Nq,Vq in enumerate(vecs):
			S[Np,Nq] = (Vp|Vq)

	evals, evecs = numpy.linalg.eigh(S)

	new_vecs = []
	for p in range(n):
		new_vec = 0
		for q in range(n):  new_vec = new_vec + evecs[q,p] * vecs[q]
		new_vecs += [ new_vec / sqrt(evals[p]) ]

	return new_vecs




class rotate_electronic_scf(object):
	def __init__(self):
		self.antisymmetrized = "antisymmetrized "	# just a convenient place to hide a printing option

	def __call__(self, h1, V2, excitations):

		# get basic information, orbital indices and their numbers
		occ_orbs = excitations._TxClasses.occ_orbs	# These lists of orbital indices ...
		vrt_orbs = excitations._TxClasses.vrt_orbs	# ... might come in uncontiguous and not numerically ordered, but only the first time.
		n_occs = len(occ_orbs)
		n_vrts = len(vrt_orbs)
		n_orbs = n_occs + n_vrts

		# read out amplitudes of all the single excitations
		singles = [ (subop._C[0],trans.index) for subop,trans in excitations ]		# looping is not recursive, so only gets singles, this feels a bit dirty (look also at hamiltonian.PT_step_conditioner and think about access functions)

		# build a list where each element is a numpy array of coefficients representing new non-ON ("dirty rotated") orbitals in the old ON basis of MOs
		coeffs = [numpy.zeros(n_orbs) for _ in range(n_orbs)]
		for p in range(n_orbs):  coeffs[p][p] = 1
		for coeff,(a,i) in singles:  coeffs[i][a] = coeff

		# build the space in which these "rotated" orbitals live, and place them "into" it
		orb_space = linear_inner_product_space(numpy_space.real_traits(n_orbs))
		orbitals = [orb_space.member(orbital) for orbital in coeffs]

		# symmetric orthonormalization of new occupied orbitals amongst themselves
		occ_orbitals = [orbitals[i] for i in occ_orbs]
		occ_orbitals = symmetric_ON(occ_orbitals)

		# Gram-Schmidt orthogonalization of new virtuals to new occupied space (which is now spanned by ON basis)
		vrt_orbitals = []
		for a in vrt_orbs:
			Va = orbitals[a]
			for Vi in occ_orbitals:
				Va = Va - Vi*(Vi|Va)	# avoid in-place modification b/c need old basis still
			vrt_orbitals += [Va]

		# symmetric orthonormalization of new virtual orbitals amongst themselves
		vrt_orbitals = symmetric_ON(vrt_orbitals)

		# build a transformation matrix from the old basis to the new basis, each column is a new MO in the old MO basis, making "0" the contraction axis
		U = numpy.zeros((n_orbs,n_orbs))
		for Ni,Vi in enumerate(occ_orbitals):
			for Np,Vp in enumerate(orbitals):
				U[Np,Ni] = (Vp|Vi)
		for Na,Va in enumerate(vrt_orbitals):
			Nx = Na + n_occs
			for Np,Vp in enumerate(orbitals):
				U[Np,Nx] = (Vp|Va)

		# Transform h1 and V2 one axis at a time
		h1_0 = h1
		h1_1 = numpy.tensordot(U, h1_0, (0,1))		#
		h1_2 = numpy.tensordot(U, h1_1, (0,1))		#
		new_h1 = h1_2					#
								#
		V2_0 = V2					# Keep in mind rules for ordering indices of the dot pdt
		V2_1 = numpy.tensordot(U, V2_0, (0,3))		#
		V2_2 = numpy.tensordot(U, V2_1, (0,3))		#
		V2_3 = numpy.tensordot(U, V2_2, (0,3))		#
		V2_4 = numpy.tensordot(U, V2_3, (0,3))		#
		new_V2 = V2_4					#

		# The rebuild above blocked by occ/vrt, rather than, say, spin.  So index lists now contiguous (and may as well be numerically ordered).
		new_fluctuations = OV_partitioning(range(0,n_occs), range(n_occs, n_orbs))

		# Compute orbital energies.   For reference   H = Sum_pq( h_pq  p* q )  +  Sum_pqrs( V_pqrs p* q* r s ),   so that   V_pqrs = (1/4) * <p(1)|<q(2)| V(1,2) [1 - P(1,2)] |s(1)>|r(2)>
		orbital_energies = []
		for p in range(n_orbs):
			energy = new_h1[p,p]
			for i in range(n_occs):			# new occs always have lowest n_occs indices now
				energy += 4 * new_V2[p,i,p,i]
			orbital_energies += [energy]

		# Compute energies for single excitations
		energy_diffs = {}
		for index in new_fluctuations.all_transition_indices().Excitations:
			a,i = index
			energy_diffs[index] = orbital_energies[a] - orbital_energies[i]

		# Finally, compute the zeroth-order and first-order, HF energy
		E0 = 0.
		for i in range(n_occs):  E0 += orbital_energies[p]
		E1 = 0.
		for i in range(n_occs):
			E1 += new_h1[i,i]
			for j in range(n_occs):			# new occs always have lowest n_occs indices now
				E1 += 2 * new_V2[i,j,i,j]

		return E0, E1, energy_diffs, new_fluctuations, new_h1, new_V2


		# I think I should be returning orbital energies and perhaps the transformation matrix ... there must be a way to define a mean-field one-body energy for a molecuel in a given state
