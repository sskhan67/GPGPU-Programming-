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
import                os
from qode.math import *
from qode      import grid_based



def Gaussian(a):
    " returns a Gaussian function with exponent a, centered at x0 "
    def g(x):  return exp(-a * x**2)
    return g



h = 2*pi
m = 1
k = 1
domain = ((-5,5),500)



S = linear_inner_product_space( grid_based.space_traits(domain) )

x = S.lin_op(   grid_based.position_operator(domain) )
n = S.lin_op( grid_based.wavenumber_operator(domain) )
p = h*n
T = p**2 / (2*m)
V = k * x**2 / 2
H = T+V

shift = 0
shifted = H + shift
guess = S.member( grid_based.function(domain,Gaussian(1)) )
guess = guess / sqrt(guess|guess)
E,v = lanczos.lowest_eigen(shifted,guess,thresh=1e-3)
print(E)
for i in range(10):
    shift = 100*(v*v)
    shifted = shifted + shift
    guess = x|v
    guess = guess / sqrt(guess|guess)
    E,v = lanczos.lowest_eigen(shifted,guess,thresh=1e-3)
    print(E)
