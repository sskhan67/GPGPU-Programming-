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
import pickle
from qode.util.parallel import resources
from qode.math          import sqrt
from qode.fermion_field import OV_partitioning
from qode.many_body.nested_operator import operator, BakerCampbellHausdorff

fluctuations = OV_partitioning([0,4],[1,2,3,5,6,7])

H_data = pickle.load(open("reference/H.pickled","rb"))
T_data = pickle.load(open("reference/T.pickled","rb"))
W_data = pickle.load(open("reference/Omega.pickled","rb"))

H     = operator.load(H_data, fluctuations)
T     = operator.load(T_data, fluctuations)
old_W = operator.load(W_data, fluctuations)

BCH = BakerCampbellHausdorff( fluctuations, resources(30) )
BCH.n_calls = 3
new_W = BCH.computeOmega(H, T, print)

print("")
print("These three numbers should match:")
print(old_W._C[0])
print(new_W._C[0])
print(open("reference/energy.txt","r").read())

norm = old_W.dot(old_W)
old_W.increment(new_W, -1)
error = sqrt(old_W.dot(old_W) / norm)

print("and the error should be small")
print("relative error =", error, "   norm =", sqrt(norm))










def apply_denominators(fock_diag, new_operator, old_operator, shift):	# Very dirty, looking at another object's privates
	if old_operator._Ex is not None:
		target = new_operator._get_Ex_block()
		for a in old_operator._Tx_indices._CrtVrt_range:
			for i in old_operator._Tx_indices._DstOcc_range:
				denom = fock_diag[a] - fock_diag[i] + shift
				target[(a,i)]._C[0] = old_operator._Ex[(a,i)]._C[0] / denom		# Dividing ... 
				apply_denominators(fock_diag, target[(a,i)], old_operator._Ex[(a,i)], denom)	# ... and then recurring misses reference, which is desired.

dT_data = pickle.load(open("reference/delT_raw.pickled","rb"))
old_dT  = operator.load(dT_data, fluctuations)

fock_diag = pickle.load(open("reference/orbital_energies.pickled","rb"))
new_dT = operator(fluctuations)
apply_denominators(fock_diag, new_dT, new_W, 0)

norm = old_dT.dot(old_dT)
old_dT.increment(new_dT, +1)			# as dumpled, old_dT already has the negative sign acted on it since we go in the opposite direction of the gradient
error = sqrt(old_dT.dot(old_dT) / norm)

print("")
print("and the error should be small here too")
print("relative error =", error, "   norm =", sqrt(norm))
