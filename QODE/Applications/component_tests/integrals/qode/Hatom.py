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
from math import *
from qode.atoms.integrals.primitive_integrals import gaussian, overlap, kinetic, charge_op

# This test digs down and uses primitive integrals

basis_fn = gaussian(1,0.28294)		# the optimal single Gaussian for an H atom.  Should give energy -0.424413
zero = (0,0,0)

S = overlap(basis_fn,zero, basis_fn,zero)

basis_fn.prefactor = 1./sqrt(S)

T = kinetic(basis_fn,zero, basis_fn,zero)
N = charge_op(basis_fn,zero, basis_fn,zero, 1,zero)

print(T+N)
print(-0.424413)
