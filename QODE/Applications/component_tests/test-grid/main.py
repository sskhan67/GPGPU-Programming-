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
from qode.math     import *
from qode          import grid_based



i      = 1j			# imaginary unit
h      = 2*pi			# Plank's constant in a.u.
m      = 1.			# mass
domain = ((-10,10),2**8)	# domain and number of grid points

def HO_centered_at(x0):		# returns a shifted potential function
    kF = 1.			# force constant
    def v(x):
        x = x - x0
        return kF * x**2 / 2.
    return v

def Gaussian(x):
    a = 1.			# exponent
    return exp(-a * x**2)



S = linear_inner_product_space( grid_based.space_traits(domain) )

V = S.lin_op( grid_based.direct_space_potential(domain,HO_centered_at(5)) )
x = S.lin_op(      grid_based.position_operator(domain) )
n = S.lin_op(    grid_based.wavenumber_operator(domain) )
p = h*n           # momentum operator (usually written hbar*k)
T = p**2 / (2*m)  # kinetic energy
H = T + V         # Hamiltonian

guess = S.member( grid_based.function(domain,Gaussian) )
guess = guess / sqrt(guess|guess)
evalue,evec = lanczos.lowest_eigen(H,guess,thresh=1e-3)
print("lowest eigenvalue = {}".format(evalue))



V = S.lin_op( grid_based.direct_space_potential(domain,HO_centered_at(0)) )

dt      = 1e-3
del_t   = 1e-2
total_t = 1e0
inner_count = int(floor(   del_t / dt ))
outer_count = int(floor( total_t / (inner_count*dt) ))

class make_frame():
    def __init__(self,domain,psi,T,V,count,dt):
        self.psi    = psi
        self.V      = V
        self.UV     = exp(-i*2*pi*V*dt)
        self.UT     = exp(-i*2*pi*T*dt)
        self.UT2    = sqrt(self.UT)
        self.count  = count
        self.t      = 0
        domain,n_points = domain
        origin,limit = domain
        dx = float(limit-origin) / n_points
        self.xvalues = [ origin+i*dx for i in range(n_points) ]
        self.figure = plot2d(ranges=(domain,(-2,18)),area=(300,200)).drawFrame().setXtics().setYtics()
        Vvalues = [ v.real for v in V.op._data ]
        self.figure.plotData(list(zip(self.xvalues,Vvalues)),markers=False,connectors=black)
    def __call__(self):
        psi = self.psi
        UV  = self.UV
        UT  = self.UT
        UT2 = self.UT2
        psi = UT2|psi
        for i in range(self.count-1):
            psi = UV|psi
            psi = UT|psi
            self.t += dt
        psi = UV|psi
        psi = UT2|psi
        self.t += dt
        self.psi = psi
        return self.plot()
    def plot(self):
        psi = self.psi
        V = self.V
        shift = (psi|V|psi).real
        xvalues = self.xvalues
        rvalues = [ shift+z.real for z in psi.v._data ]
        ivalues = [ shift+z.imag for z in psi.v._data ]
        pvalues = [ shift+abs(z) for z in psi.v._data ]
        nvalues = [ shift-abs(z) for z in psi.v._data ]
        figure = self.figure()
        figure.plotData(list(zip(xvalues,pvalues)),markers=False,connectors=black)
        figure.plotData(list(zip(xvalues,nvalues)),markers=False,connectors=black)
        figure.plotData(list(zip(xvalues,rvalues)),markers=False,connectors=red)
        figure.plotData(list(zip(xvalues,ivalues)),markers=False,connectors=green)
        return figure
frame = make_frame(domain,evec,T,V,inner_count,dt)

animate.svg("psi",frame,outer_count,10)
