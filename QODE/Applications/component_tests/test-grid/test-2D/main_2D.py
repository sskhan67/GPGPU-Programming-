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
import numpy


i      = 1j			                # imaginary unit
h      = 2*pi		                	# Plank's constant in a.u.
m      = 1.             		       	# mass
hbar   = h/(2*pi)
x_points = 2*2**8
y_points = 2**8
domain = ( ( (-20,20), (-10,10) ), ((x_points), (y_points) ) )	# domain and number of grid points
kF     = 1.


def Gaussian(x):
    a = 1.			# exponent
    return exp(-a * (x-5)**2)


def Gaussian_2D(x,y):
   alpha_x = ( kF*m/(hbar**2) )** 0.5
   alpha_y = ( kF*m/(hbar**2) )** 0.5
   Norm_x  = (alpha_x/pi)**0.25
   Norm_y  = (alpha_y/pi)**0.25
   expon_x = exp( (-1/2)*alpha_x* (x-5)**2) 
   expon_y = exp( (-1/2)*alpha_y* (y-5)**2) 
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

S = linear_inner_product_space( grid_based.space_traits(domain) )

# x = S.lin_op(      grid_based.position_operator(domain) )
# n = S.lin_op(    grid_based.wavenumber_operator(domain) )
# p = h*n           # momentum operator (usually written hbar*k)
# T = p**2 / (2*m)  # kinetic energy

V = S.lin_op( grid_based.direct_space_potential(domain,HO_2D_centered_at(5,5)) )
L = S.lin_op( grid_based.laplacian_operator(domain) )
T = -(hbar**2)/(2*m) * L   #kinetic energy operator 
H = T + V                  # Hamiltonian

guess = S.member( grid_based.function(domain,Gaussian_2D) )

guess = guess / sqrt(guess|guess)
V_average = (guess|V|guess)
T_average = (guess|T|guess)
Total_energy = V_average + T_average
print("Kinetic Energy",   T_average)
print("Potential Energy", V_average)
print("Total Energy", Total_energy)

# evalue,evec = lanczos.lowest_eigen(H,guess,thresh=1e-3)
# print("lowest eigenvalue = {}".format(evalue))

V = S.lin_op( grid_based.direct_space_potential(domain,HO_2D_centered_at(0,0)) )

dt      = 1e-3
del_t   = 1e-2
total_t = 1e0
inner_count = int(floor(   del_t / dt ))
outer_count = int(floor( total_t / (inner_count*dt) ))

# This sets are for testing purposes only, use above values for FULL animation
outer_count = 4                 
inner_count = 5

def create_gnuplot_file(xvalues, yvalues,Vvalue):
        output_str  = ""
        output_str  = "{} ".format(xvalues)
        output_str += "{} ".format(yvalues)
        output_str += "{}\n".format(Vvalue)
        return(output_str)

def spaces():
        output_str = ""
        output_str += "\n"
        output_str += "\n"
        output_str += "## NEW FILE ###\n"
        return(output_str)


class make_frame():
    def __init__(self,domain,psi,T,V,count,dt):
        self.psi    = psi
        self.V      = V
        self.UV     = exp(-i*2*pi*V*dt)
        self.UT     = exp(-i*2*pi*T*dt)
        self.UT2    = sqrt(self.UT)
        self.count  = count
        self.t      = 0
        self.index  = 0
        domain,n_points = domain
        self.dimension = len(domain)
        if self.dimension == 1:
            origin,limit = domain
            dx = float(limit-origin) / n_points
            self.xvalues = [ origin+i*dx for i in range(n_points) ]
            self.figure = plot2d(ranges=(domain,(-2,18)),area=(300,200)).drawFrame().setXtics().setYtics()
            Vvalues = [ v.real for v in V.op._data ]
            self.figure.plotData(list(zip(self.xvalues,Vvalues)),markers=False,connectors=green)   #plot y_values??
        if self.dimension == 2:
            x_domain,y_domain = domain
            x_origin,x_limit = x_domain
            y_origin,y_limit = y_domain
            xn_points, yn_points = n_points
            dx = float(x_limit-x_origin) / xn_points
            dy = float(y_limit-y_origin) / yn_points
            self.xvalues = [ x_origin+i*dx for i in range(xn_points) ]
            self.yvalues = [ y_origin+i*dy for i in range(yn_points) ]
            Vvalues = V.op._data.real.tolist()
            f = open('V_data.txt', 'a')
            dimx = len(Vvalues)
            dimy = len(Vvalues[0])
            for row in range(dimx):
             for col in range(dimy):
              f.write( create_gnuplot_file(self.xvalues[row], self.yvalues[col], Vvalues[row][col]) )
            f.close()
            # data2bmp('Vvalues',Vvalues, Zscale=0.008)   # creates a graph of potential used
            self.figure = composite()
            self.x_origin = x_origin
            self.x_limit  = x_limit
            self.y_origin = y_origin
            self.y_limit  = y_limit
            #self.figure += image(filename='Vvalues.bmp',translate=(100,-100))
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
#       shift = (psi|V|psi).real
        shift = 0
        if self.dimension == 1:
            xvalues = self.xvalues
            rvalues = [ shift+z.real for z in psi.v._data ]
            ivalues = [ shift+z.imag for z in psi.v._data ]
            pvalues = [ shift+abs(z) for z in psi.v._data ]
            nvalues = [ shift-abs(z) for z in psi.v._data ]
            figure.plotData(list(zip(xvalues,pvalues)),markers=False,connectors=black)
            figure.plotData(list(zip(xvalues,nvalues)),markers=False,connectors=black)
            figure.plotData(list(zip(xvalues,rvalues)),markers=False,connectors=red)
            figure.plotData(list(zip(yvalues,rvalues)),markers=False,connectors=red)
            figure.plotData(list(zip(xvalues,ivalues)),markers=False,connectors=green)
            figure.plotData(list(zip(yvalues,ivalues)),markers=False,connectors=green)
        if self.dimension == 2:
            xvalues = self.xvalues
            yvalues = self.yvalues
            rvalues = psi.v._data.real.tolist()
            #  The next lines create files of each time frame of wavefunction
            #  f = open('r_data.txt', 'a')
            #  rdimx = len(rvalues)
            #  rdimy = len(rvalues[0])
            #  print(rdimx, 'rdimx')
            #  print(rdimy, 'rdimy')
            #  for row in range(rdimx):
            #   for col in range(rdimy):
            #    f.write( create_gnuplot_file(self.xvalues[row], self.yvalues[col], rvalues[row][col]) )
            #  f.write( spaces() )
            #  f.close()
            filestem = 'values/values_{}'.format(str(self.index).zfill(3))
            # filestem = 'kvalues/kvalues_{}'.format(str(self.index).zfill(3))
            rotated_data = numpy.rot90( psi.v._data )  # modified line
            # kvalues = numpy.fft.fft2(psi.v._data)
            data2bmp(filestem, rotated_data, Zscale=1, colormap=complex2color(0.1) )
            # data2bmp(filestem, kvalues, Zscale=1, colormap=complex2color(0.1) )
            #self.figure += image(filename='rvalues.bmp', translate=(100,-100))
            #self.figure += image(filename='ivalues.bmp', translate=(100,-100))
            figure = self.figure()
            figure += image(filename=filestem+'.bmp',scale=3)
            aspect = 300*((self.y_limit - self.y_origin)/(self.x_limit - self.x_origin)) 
            coords = plot2d(ranges=((self.x_origin,self.x_limit),(self.y_origin,self.y_limit)),area=(300,aspect)).drawFrame()
            coords.setXtics(connectors=white,textform=text("{:.0f}",translate=(0,-10)))
            coords.setYtics(connectors=white,textform=text("{:.0f}",translate=(-10,0)))
            figure += coords(translate=(0,aspect))
            figure += text('$t$={:10.3f}'.format(self.t),color=white,  translate=(20,aspect-10),alignment=(left,center), typeface='latex-ComputerModern' )
        figure = figure(translate=(100,-600))
        self.index += 1
        return figure
frame = make_frame(domain,guess,T,V,inner_count,dt)

animate.svg("psi_2D",frame,outer_count,10)
