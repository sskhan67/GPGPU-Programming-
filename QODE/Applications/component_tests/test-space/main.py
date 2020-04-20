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
from qode      import util
from qode.math import linear_inner_product_space as space
from qode.math import numpy_space

S = space( numpy_space.real_traits(10) )

v = util.basis_vector(2,10)
w = util.basis_vector(3,10)

a = S.member(v)
b = S.member(w)

print( (a|b) )
print( (a|a) )
print( (b|b) )
