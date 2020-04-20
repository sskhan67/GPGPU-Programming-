#    (C) Copyright 2018 Yuhong Liu
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
import pickle

np.set_printoptions(linewidth=270,threshold=np.nan,precision=2)


def spin_to_spat(index, n_spat):
	if index < n_spat:
		return index
	else:
		return index - n_spat

def same_spin(orb1, orb2, n_spat):
	if orb1 < n_spat:
		if orb2 < n_spat:
			return True
		else:
			return False
	else:
		if orb2 < n_spat:
			return False
		else:
			return True

def H_element(bra, hmat, Vmat, ket):
	n_spat = hmat.shape[0]
	print("num spatial  =", n_spat)
	if bra == ket:
		# print('same')
		res = 0.0
		for m in ket:
			print('adding h[%d] = %f' %(m, hmat[spin_to_spat(m,n_spat),spin_to_spat(m,n_spat)]))
			res += hmat[spin_to_spat(m,n_spat),spin_to_spat(m,n_spat)]
		for m in ket:
			for n in ket:
				print('adding  V[%d,%d,%d,%d] = %f' %(m,n,m,n,  Vmat[spin_to_spat(m,n_spat),spin_to_spat(n,n_spat),spin_to_spat(m,n_spat),spin_to_spat(n,n_spat)]))
				res += 0.5 * Vmat[spin_to_spat(m,n_spat),spin_to_spat(n,n_spat),spin_to_spat(m,n_spat),spin_to_spat(n,n_spat)]
				if same_spin(m,n,n_spat):
					print('adding -V[%d,%d,%d,%d] = %f' %(m,n,n,m, -Vmat[spin_to_spat(m,n_spat),spin_to_spat(n,n_spat),spin_to_spat(n,n_spat),spin_to_spat(m,n_spat)]))
					res += -0.5 * Vmat[spin_to_spat(m,n_spat),spin_to_spat(n,n_spat),spin_to_spat(n,n_spat),spin_to_spat(m,n_spat)]
		return res
	else:
		print('not same')
		sbra = set(bra)
		sket = set(ket)
		n_diff = len( list(sbra ^ sket) )  // 2
		if n_diff == 1:
			...
		if n_diff == 2:
			...


if __name__ == '__main__':
	# Be dimer/ HF energy/ 6-31G
	# Psi4 energy           = -29.135584149665647
	# Boris's script energy = -29.1244915865
	# Energy matches Psi4 False
	# Nuclear repulsion energy = 3.2564751297846146

	# H = pickle.load( open('H.p','rb') )
	# print(H)

	hmat = np.load('biorthogonal_core_matrix_dimer.npy')
	Vmat = np.load('biorthogonal_twobody_matrix_dimer.npy')

	# print(hmat[3,3]+hmat[6,6]+hmat[9,9]+hmat[12,12])

	# bra = [0,3,9,12,18,21,27,30]
	# ket = [0,3,9,12,18,21,27,30]
	bra = [0,1,9,10,18,19,27,28]
	ket = [0,1,9,10,18,19,27,28]

	# print(hmat)

	value = H_element(bra, hmat, Vmat, ket)
	print('my H element =',value)


	# H_element(bra, hmat, Vmat, [0,9,10,12])
	# H_element(bra, hmat, Vmat, [0,9,10,11])
