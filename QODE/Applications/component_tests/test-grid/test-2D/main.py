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
from potential     import HO_centered_at
from potential     import barrier_at

i      = 1j			# imaginary unit
h      = 2*pi			# Plank's constant in a.u.
m      = 1.			# mass
hbar   = h/(2*pi)               
domain = ((-10,10), (2**8) )	# domain and number of grid points

def Gaussian(x):
    a = 1.			# exponent
    return exp(-a * (x)**2)

def Gaussian_1D(x):             # Analytical solution for ground state of HO
   kF = 1
   alpha          = ( kF*m/(hbar**2) )** 0.5 
   Norm_constant  = (alpha/pi)**0.25
   wavefunc       = exp( (-1/2)*alpha* (x-5)**2)
   return  Norm_constant * wavefunc   


def wavepacket(x):
    v = 3
    a = 1
    return exp( -a*(x+5)**2)*exp(i*pi*v*x)

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

S = linear_inner_product_space( grid_based.space_traits(domain) )

V = S.lin_op( grid_based.direct_space_potential(domain,HO_centered_at(5)) )
x = S.lin_op(      grid_based.position_operator(domain) )
n = S.lin_op(    grid_based.wavenumber_operator(domain) )
p = h*n           # momentum operator (usually written hbar*k)
T = p**2 / (2*m)  # kinetic energy
H = T + V         # Hamiltonian

guess = S.member( grid_based.function(domain,Gaussian_1D) )
guess = guess / sqrt(guess|guess)

# OLD SYNTAX
# evalue, evec = lanczos.lowest_eigen(H,guess,thresh=1e-3)


# NEW SINTAX
evpairs = list( lanczos.lowest_eigen(H,[guess],thresh=1e-3) )
evalue = evpairs[0][0]
evec   = evpairs[0][1].v
print(evec)
print(guess)
print("lowest eigenvalue = {}".format(evalue))

# V_average = ( evec|V|evec )
# T_average = ( evec|T|evec )
V_anal    = (guess|V|guess )
T_anal    = (guess|T|guess )
# Total_energy      = V_average + T_average
Total_energy_anal = V_anal + T_anal
# print("Kinetic Energy",        T_average         )
# print("Potential Energy",      V_average         )
# print("Total Energy",          Total_energy      )
print("Kinetic analytic",      T_anal            )
print("Potential analytic",    V_anal            )
print("Total energy analytic", Total_energy_anal )


V = S.lin_op( grid_based.direct_space_potential(domain,HO_centered_at(0)) )

dt      = 1e-3
del_t   = 1e-2
total_t = 1e0
inner_count = int(floor(   del_t / dt ))
outer_count = int(floor( total_t / (inner_count*dt) ))

# This arbitrary number is for testing purposes only
outer_count = 50



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
frame = make_frame(domain,guess,T,V,inner_count,dt)

animate.svg("psi_1D",frame,outer_count,10)
