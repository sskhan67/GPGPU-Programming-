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
from pytoon        import *
from pytoon.macros import *
from pytoon.raster import *
from qode.math     import *
from qode          import grid_based


# this file contains all paramters for this directory

i      = 1j                                     # imaginary unit
h      = 2*pi                                   # Plank's constant in a.u.
m      = 1.                                     # mass
hbar   = h/(2*pi)
domain = ( ( (-40,40), (-5,5) ), ( ( 4*(2**7) ), ( (2**6) ) ) )     # domain and number of grid points
kF     = 7.
m1     = 1.
m2     = 2.
radius = 3.
height = 3.
