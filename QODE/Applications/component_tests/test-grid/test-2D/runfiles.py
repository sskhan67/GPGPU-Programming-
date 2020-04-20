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
from potential     import *
from complex2color import *
from wavefunction  import * 
from parameters    import *
import numpy
import os


# x = S.lin_op(      grid_based.position_operator(domain) )
# n = S.lin_op(    grid_based.wavenumber_operator(domain) )
# p = h*n           # momentum operator (usually written hbar*k)
# T = p**2 / (2*m)  # kinetic energy
# L = S.lin_op( grid_based.laplacian_operator(domain) )
# T = -(hbar**2)/(2*m) * L   #kinetic energy operator 
# L1 = S.lin_op( grid_based.laplacian_operator_x1(domain) )
# L2 = S.lin_op( grid_based.laplacian_operator_x2(domain) )
# T1 = -(hbar**2)/(2*m1) * L1   #kinetic energy operator particle 1
# T2 = -(hbar**2)/(2*m2) * L2   #kinetic energy operator particle 2
# T = T2 + T1



#V_average = (guess|V|guess)
#T_average = (guess|T|guess)
#Total_energy = V_average + T_average
#print("Kinetic Energy",   T_average)
#print("Potential Energy", V_average)
#print("Total Energy", Total_energy)

# evalue,evec = lanczos.lowest_eigen(H,guess,thresh=1e-3)
# print("lowest eigenvalue = {}".format(evalue))



# V = S.lin_op( grid_based.direct_space_potential(domain,two_particle_barrier_HO(height,radius,0)) )


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
            Vvalues = numpy.rot90( V.op._data ).real.tolist()
            f = open('V_data.txt', 'a')
            dim = len(Vvalues)
            for row in range(dim):
             for col in range(dim):
              f.write( create_gnuplot_file(self.xvalues[row], self.yvalues[col], Vvalues[row][col]) )
            f.close()
            data2bmp('Vvalues',Vvalues, Zscale=0.008)
            self.figure = composite()
            #self.figure += image(filename='Vvalues.bmp',translate=(100,-100))
    def __call__(self):
        fig = self.plot()
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
        return fig
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
            #  ivalues = psi.v._data.imag.tolist() 
            #  pvalues =    numpy.abs( psi.v._data ).tolist() 
            #  nvalues = ( -numpy.abs( psi.v._data ) ).tolist()
            #  data2bmp('rvalues', rvalues, Zscale=0.008, colormap=complex2color_log() )
            #  data2bmp('ivalues', ivalues, Zscale=0.008, colormap=complex2color_log() )
            # f = open('r_data.txt', 'a')
            # rdimx = len(rvalues)
            # rdimy = len(rvalues[0])
            # for row in range(rdimx):
            #   for col in range(rdimy):
            #    f.write( create_gnuplot_file(self.xvalues[row], self.yvalues[col], rvalues[row][col]) )
            # f.write( spaces() )
            # f.close()
            filestem = 'kvalues/kvalues_{}'.format(str(self.index).zfill(3))
            # filestem  = 'values/values_{}'.format(str(self.index).zfill(3))
            kvalues = numpy.fft.fft2(   psi.v._data )
            rotated_data = numpy.rot90( kvalues )
            data2bmp(filestem, rotated_data, Zscale=1, colormap=complex2color(0.1))
            # data2bmp(filestem,  rotated_data, Zscale=1, colormap=complex2color(0.1))
            # self.figure += image(filename='rvalues.bmp', translate=(100,-100))
            # self.figure += image(filename='ivalues.bmp', translate=(100,-100))
            figure = self.figure()
            # figure += image(filename=filestem+'.bmp',scale=3)
            figure += image(filename=filestem+'.bmp',scale=3) 
            coords = plot2d(ranges=((-40,40),(-10,10)),area=(300,75)).drawFrame()
            coords.setXtics(connectors=white,textform=text("{:.0f}",xheight=10, translate=(0,-10)))
            coords.setYtics(connectors=white,textform=text("{:.0f}",xheight=10, translate=(-10,0)))
            figure += coords(translate=(0,75))
            figure += line(begin=(0,37.5), end=(300,37.5), color=white, weight=0.25 )
            figure += text('$t$={:10.3f}'.format(self.t),color=white, xheight=10, translate=(10,65),alignment=(left,center), typeface='latex-ComputerModern' )
        figure = figure(translate=(100,-100))
        self.index += 1
        return figure

def setup(m2,height,kF):
	S = linear_inner_product_space( grid_based.space_traits(domain) )
	V = S.lin_op( grid_based.direct_space_potential(domain,two_particle_barrier_HO(height,radius,0,kF) ) )
	L = S.lin_op( grid_based.laplacian_operator_x1_x2(domain, m1, m2) )
	T = -(hbar**2)/2 * L
	H = T + V                  # Hamiltonian
	guess = S.member( grid_based.function(domain,composite_wfn) )
	guess = guess / sqrt(guess|guess)
	frame = make_frame(domain,guess,T,V,inner_count,dt)
	animate.svg("psi",frame,outer_count,10)
	output_str  = ''
	output_str += 'trail_kF_{}'.format(kF)
	output_str += '__m2_{}'.format(m2)
	output_str += '__h_{}'.format(height)
	os.system('mkdir ' + 'two_particles_animations/' + output_str ) 
	os.system('cp psi.svg ' + 'two_particles_animations/' + output_str + '/')
	os.system('cp Vvalues.bmp ' + 'two_particles_animations/' + output_str + '/')
	os.system('cp -r  values/ two_particles_animations/' + output_str + '/')      
	os.system('cp -r kvalues/ two_particles_animations/' + output_str + '/')      

dt      = 1e-3
del_t   = 1e-2
total_t = 1e0
inner_count = int(floor(   del_t / dt ))
outer_count = int(floor( total_t / (inner_count*dt) ))
outer_count = 5 

for kF in [5,10]:
    for m2 in [1,4]:
         for height in [10, 5]:
             setup(m2, height, kF)
