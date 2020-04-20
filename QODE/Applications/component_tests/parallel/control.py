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
import time
import numpy
import qode



class op_vec(object):
	def __init__(self, out, in_vec):
		self.out    = out
		self.in_vec = in_vec
	def __call__(self,mat):
		self.out += mat.dot(self.in_vec)

in_vec = qode.util.random_matrix(1000,seed=2001).dot(qode.util.basis_vector(7,1000))
mat    = qode.util.random_matrix(1000,seed=1984)



start = time.time()

out = qode.util.zero_vector(1000)
f = op_vec(out, in_vec)

for _ in range(100):  f(mat)

end = time.time()
runtime = end - start

print("control imported")
