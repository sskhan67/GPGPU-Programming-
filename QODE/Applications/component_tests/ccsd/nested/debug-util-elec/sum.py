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

op_files = [\
"from_EED_E.pickled",
"from_EED_EE.pickled",
"from_EFF_E.pickled",
"from_EFF_EE.pickled",
"from_ED_E.pickled",
"from_ED_EE.pickled",
"from_DD_E.pickled",
"from_DD_EE.pickled",
"from_EF_E.pickled",
"from_EF_EE.pickled",
"from_F_E.pickled",
"from_F_EE.pickled",
"from_D_E.pickled",
"from_D_EE.pickled",
"from_FD_E.pickled",
"from_FD_EE.pickled",
"from_EFD_E.pickled",
"from_EFD_EE.pickled",
"from_FF_E.pickled",
"from_FF_EE.pickled"
]



order = int(sys.argv[1])

fluctuations = OV_partitioning([0,4],[1,2,3,5,6,7])

op_sum = operator(fluctuations)

for op_file in op_files:
	op_path = "masked/X_{}/{}".format(order,op_file)
	op_data = pickle.load(open(op_path,"rb"))
	op = operator.load(op_data, fluctuations)
	op_sum.increment(op)

pickle.dump(op_sum.dump(),open("masked/X_{}/sum.pickled".format(order),"wb"))
