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
from qode.math          import sqrt
from qode.fermion_field import OV_partitioning
from qode.many_body.nested_operator import operator, BakerCampbellHausdorff, mask
from the_masks import all_masks



fluctuations = OV_partitioning([0,4],[1,2,3,5,6,7])

op_file = sys.argv[1]

op_data = pickle.load(open(op_file,"rb"))

if len(sys.argv)==3:
	the_mask = all_masks[sys.argv[2]]
	op_data = the_mask(operator.load(op_data,fluctuations)).dump()

for indices,amplitude in op_data:
	print("{: .2e}".format(amplitude),end="   ")
	for index in indices:
		itype,(p,q) = index
		print("{}({},{})".format(itype,p,q), end=" ")
	print()
