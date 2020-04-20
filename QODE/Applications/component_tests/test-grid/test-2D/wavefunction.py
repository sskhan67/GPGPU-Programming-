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
from potential     import HO_2D_centered_at
from potential     import barrier_2D_at
from complex2color import *
from parameters    import *
import numpy



def Gaussian(x):
    a = 1.                      # exponent
    return exp(-a * (x-5)**2)


def Gaussian_2D(x,y):
   alpha_x = 10* ( kF*m/(hbar**2) )** 0.5
   alpha_y = ( kF*m/(hbar**2) )** 0.5
   Norm_x  = (alpha_x/pi)**0.25
   Norm_y  = (alpha_y/pi)**0.25
   expon_x = exp( (-1/2)*alpha_x* (x)**2)
   expon_y = exp( (-1/2)*alpha_y* (y)**2)
   return Norm_x*expon_x * Norm_y*expon_y

   # a = ((2)**0.5)/4
   # sigma_squared = hbar/ (2 * (kF*m)**0.5 )
   # a = 1/( 4 * sigma_squared )
   # return exp( -a*( (x)**2 + (y)**2 ) )

def wavepacket(x):
    v = 3
    a = 1
    return exp( -a*(x+5)**2)*exp(i*pi*v*x)

def wavepacket_2D(x,y):
    v = 3
    a = 1
    return exp( -a*( (x+5)**2 + (y)**2 ) )*exp(i*pi*v*x)

def wavepacket_lin_comb(x):
    v = 3                       # wavenumber
    a = 1                       # exponent
    for index in range(10,31):
        k        = index/10
        wavPkt   = 0
        gausCoef = exp(-a * ( k-2 )**2 )
        planeW   = exp(i * k * (x+5) )
        wavPkt  += planeW
    return gausCoef*wavPkt

def composite_wfn(kF,mass2):
    a = 1/100.                      # exponent
    v = 1.5    #1.5 previously
    alpha = ( kF*mass2/( hbar**2 ) )**0.5
    def wfn(x1,x2):
        x1 = x1+20
        gaussian    = exp( (-1/2)* alpha * (x2)**2)
        wavepacket  = exp( -a*(x1)**2)*exp(i*pi*v*x1)
        wavepacket2 = exp( -a*(x1-120)**2)*exp(i*pi*v*x1)
        return (gaussian * wavepacket) + (gaussian * wavepacket2)
    return wfn
