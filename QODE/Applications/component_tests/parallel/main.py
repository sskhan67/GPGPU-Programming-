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
import math
import qode
import control
from control import mat, in_vec

# The speedup here is actually <1, but that is just because the inputs we are sending are about the same size as the computation.
# Originally used a bunch of locally generated random matrices, but concurrent calls to random slowed things way down.
# Would be more fun/realistic to wrap func in a class so that it stores its own instance of random.Random, then just send seeds.


def func(i,mat,vec):
	print(i)
	return mat.dot(vec)

def value_adder(out,inc):
	out += inc
	return out

start = time.time()

A = qode.util.parallel.aggregate(func, (in_vec,), value_adder, 10)
in_mat = [(i,mat) for i in range(100)]
out = qode.util.zero_vector(1000)
out = A.run(out,in_mat)

end = time.time()
runtime = end - start

print("experiment done")



print("control runtime =", control.runtime)
print("speedup =", control.runtime/runtime)
print()

true_ctrl = 100 * mat.dot(in_vec)
ctrl = control.out
expt = out

print("vector norm =", math.sqrt(true_ctrl.dot(true_ctrl)))
diff = ctrl - true_ctrl
error = math.sqrt( diff.dot(diff) / ctrl.dot(ctrl) )
print("control relative error =", error)
diff = expt - true_ctrl
error = math.sqrt( diff.dot(diff) / ctrl.dot(ctrl) )
print("experiment relative error =", error)
