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
from qode.many_body.nested_operator import operator, BakerCampbellHausdorff

fluctuations = OV_partitioning([0,4],[1,2,3,5,6,7])

op1_file = sys.argv[1]
op2_file = sys.argv[2]

op1_data = pickle.load(open(op1_file,"rb"))
op2_data = pickle.load(open(op2_file,"rb"))

op1 = operator.load(op1_data, fluctuations)
op2 = operator.load(op2_data, fluctuations)

norm = op1.dot(op1)
if norm>0:
	op1.increment(op2, -1)
	error = sqrt(op1.dot(op1) / norm)
	print("relative error =", error, "   norm =", sqrt(norm))
else:
	norm = op2.dot(op2)
	print("vanishing norm for one operator.  norm of other=", sqrt(norm))
