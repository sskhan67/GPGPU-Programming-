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
import sys
import pickle
from qode.util.parallel import resources
from qode.fermion_field import OV_partitioning
from qode.many_body.nested_operator import operator, mask, commute
from qode.many_body.nested_operator.baker_campbell_hausdorff.BCHtruncation import BCHtruncation
from the_masks import Xmasks, Tmasks



fluctuations = OV_partitioning([0,4],[1,2,3,5,6,7])

order = int(sys.argv[1])

wisdom = BCHtruncation(2,2).masks(order)

X_data = pickle.load(open("masked/X_{}/total.pickled".format(order-1),"rb"))
T_data = pickle.load(open("masked/T.pickled","rb"))

X = operator.load(X_data, fluctuations)
T = operator.load(T_data, fluctuations)

commutator = commute( X, T, resources(30) )
commutator.scale(1./order)

result = wisdom.to_keep(commutator).copy()
result.increment(wisdom.excitations(commutator))

pickle.dump(result.dump(), open("masked/X_{}/check.pickled".format(order),"wb"))

for Xlabel,Xmask in Xmasks.items():
	for Tlabel,Tmask in Tmasks.items():
		commutator = commute( Xmask(X), Tmask(T), resources(30) )
		commutator.scale(1./order)
		result = wisdom.to_keep(commutator).copy()
		result.increment(wisdom.excitations(commutator))
		pickle.dump(result.dump(), open("masked/X_{}/from_{}_{}.pickled".format(order,Xlabel,Tlabel),"wb"))
