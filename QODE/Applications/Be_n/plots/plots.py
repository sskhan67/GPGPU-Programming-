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
import numpy
from pytoon        import *
from pytoon.macros import *

def to_mEh(nrg):  return 1000*nrg
def to_min(t):    return t/60.

tiger = "#F47920"

purple_cross = composite()
purple_cross += line(begin=(-0.5,-0.5),end=(0.5,0.5),color=purple,weight=1/3.)
purple_cross += line(begin=(0.5,-0.5),end=(-0.5,0.5),color=purple,weight=1/3.)
purple_cross = purple_cross(scale=3)
blue_cross = composite()
blue_cross += line(begin=(-0.5,-0.5),end=(0.5,0.5),color=blue,weight=1/3.)
blue_cross += line(begin=(0.5,-0.5),end=(-0.5,0.5),color=blue,weight=1/3.)
blue_cross = blue_cross(scale=3)
tiger_cross = composite()
tiger_cross += line(begin=(-0.5,-0.5),end=(0.5,0.5),color=tiger,weight=1/3.)
tiger_cross += line(begin=(0.5,-0.5),end=(-0.5,0.5),color=tiger,weight=1/3.)
tiger_cross = tiger_cross(scale=3)
red_cross = composite()
red_cross += line(begin=(-0.5,-0.5),end=(0.5,0.5),color=red,weight=1/3.)
red_cross += line(begin=(0.5,-0.5),end=(-0.5,0.5),color=red,weight=1/3.)
red_cross = red_cross(scale=3)
tiger_circle = circle(radius=3*0.5,color=tiger,fill=None,weight=1)



from data import be2_Xccsd
from data import be2_Xccsd_noCT
from data import be2_ccsd
from data import be2_ccsd_t
from data import be2_fci

from data import be3_Xccsd
from data import be3_Xfci
from data import be3_Xfci_7
from data import be3_ccsd
from data import be3_ccsd_t
from data import be3_fci

from data import beN_Xccsd_timing
from data import beN_ccsd_timing
from data import beN_t_timing

from data import beN_Xccsd_welldepth
from data import beN_ccsd_welldepth
from data import beN_ccsd_t_welldepth



textX = text(xheight=15)


dimer = plot2d(ranges=((2.99,9.99),(-7e-4,1e-4)), area=(300,200)).drawFrame()
dimer.plotData(be2_fci.data,    connectors=(black,3),markers=False)
dimer.plotData(be2_ccsd_t.data, connectors=line(color=green,style=dotted),    markers=False)
dimer.plotData(be2_ccsd.data,   connectors=line(color=blue,style=dashed),     markers=False)
dimer.plotData(be2_Xccsd.data,  connectors=False,    markers=tiger_cross(scale=1.5))
dimer.plotData(be2_Xccsd_noCT.data,  connectors=False,    markers=tiger_circle(scale=1.5))
dimer.setXtics(periodicity=1,   textform=textX("{:.0f}",alignment=(center,top),  translate=(0,-10)))
dimer.setYtics(periodicity=2e-4,textform=(textX("{:.1f}",alignment=(right,center),translate=(-10,-0)),to_mEh))
dimer += textX("\\textbf{Be}$_\\textbf{2}$", translate=(150,210),alignment=(center,bottom))
#dimer += textX("$E_\\text{h}$", translate=(-95,100),alignment=(right,center))
#dimer += textX("\\AA", translate=(150,-35),alignment=(center,top))
#dimer.pdf("be2",texlabels="texlabels/2/text*.svg")

trimer = plot2d(ranges=((3.01,10.01),(-7e-4,1e-4)), area=(300,200)).drawFrame()
trimer.plotData(be3_fci.data,       connectors=(black,3),markers=False)
trimer.plotData(be3_ccsd_t.data,    connectors=line(color=green,style=dotted),    markers=False)
trimer.plotData(be3_ccsd.data,      connectors=line(color=blue,style=dashed),     markers=False)
trimer.plotData(be3_Xfci.data,      connectors=False,    markers=tiger_circle(scale=1.5))
trimer.plotData(be3_Xccsd.data,     connectors=False,    markers=tiger_cross(scale=1.5))
#trimer.plotData(be2_fci.inc_trimer, connectors=yellow,markers=False)
#trimer.plotData(be3_Xfci_7.data,    connectors=False,    markers=tiger_cross)
trimer.setXtics(periodicity=1,   textform=textX("{:.0f}",alignment=(center,top),  translate=(0,-10)))
#trimer.setYtics(periodicity=2e-4,textform=(textX("{:.1f}",alignment=(right,center),translate=(-10,-0)),to_mEh))
trimer.setYtics(periodicity=2e-4,textform=(textX("",alignment=(right,center),translate=(-10,-0)),to_mEh))
trimer += textX("\\textbf{Be}$_\\textbf{3}$", translate=(150,210),alignment=(center,bottom))
#trimer += textX("$E_\\text{h}$", translate=(-95,100),alignment=(right,center))
#trimer += textX("\\AA", translate=(150,-35),alignment=(center,top))
#trimer.pdf("be3",texlabels="texlabels/3/text*.svg")

figure1 = composite()
figure1 += polygon(points=[(-80,30),(635,30),(635,-255),(-80,-255)],fill=white,weight=0)
figure1 += dimer()
figure1 += trimer(translate=(320,0))
figure1 += textX("Internuclear separation / \\AA", translate=(310,-230),alignment=(center,top))
figure1 += textX("Bond potential / m\\textit{E}$_\\text{h}$", translate=(-65,-100),alignment=(center,center),rotate=90)
w,h,x,y,o = 140,90,200,85,4
figure1 += polygon(points=[(0,0),(0,-h),(w,-h),(w,0)],translate=(x+o,-(y+o)),fill=gray,weight=0)
figure1 += polygon(points=[(0,0),(0,-h),(w,-h),(w,0)],translate=(x,-y),fill=white)
x0,x1,y0,dx,dy = x+10,x+40,y+15,7,20
figure1 += line(color=blue,style=dashed,begin=(x0,-y0),end=(x1,-y0))
figure1 += textX("CCSD",translate=(x1+5,-y0),alignment=(left,center))
figure1 += line(color=green,style=dotted,begin=(x0,-(y0+dy)),end=(x1,-(y0+dy)))
figure1 += textX("CCSD(T)",translate=(x1+5,-(y0+dy)),alignment=(left,center))
figure1 += line(weight=3,begin=(x0,-(y0+2*dy)),end=(x1,-(y0+2*dy)))
figure1 += textX("FCI",translate=(x1+5,-(y0+2*dy)),alignment=(left,center))
figure1 += tiger_cross(translate=(x0+0.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += tiger_cross(translate=(x0+1.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += tiger_cross(translate=(x0+2.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += tiger_cross(translate=(x0+3.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += textX("X2-CCSD",translate=(x1+5,-(y0+3*dy)),alignment=(left,center))
x0 = x0 - 190
x1 = x1 - 190
y0 = y0 - 8
figure1 += tiger_circle(translate=(x0+0.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += tiger_circle(translate=(x0+1.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += tiger_circle(translate=(x0+2.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += tiger_circle(translate=(x0+3.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += textX("X2-CCSD without \\\\ charge transfer",translate=(x1+5,-(y0+3*dy)),alignment=(left,center))
x0 = x0 + 425
x1 = x1 + 425
figure1 += tiger_circle(translate=(x0+0.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += tiger_circle(translate=(x0+1.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += tiger_circle(translate=(x0+2.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += tiger_circle(translate=(x0+3.5*dx,-(y0+3*dy)),scale=1.5)
figure1 += textX("X2-FCI \\\\ (= X2-CCSDT)",translate=(x1+5,-(y0+3*dy)),alignment=(left,center))
figure1.pdf("PEcurves",texlabels="texlabels/PE/text*.svg")














fit_data_x = [log(x) for x,y in beN_Xccsd_timing.data[-10:]]
fit_data_y = [log(y) for x,y in beN_Xccsd_timing.data[-10:]]
fit = numpy.polyfit(fit_data_x,fit_data_y,1)
XccsdFit = exp(fit[1]), fit[0]

fit_data_x = [log(x) for x,y in beN_ccsd_timing.data[-10:]]
fit_data_y = [log(y) for x,y in beN_ccsd_timing.data[-10:]]
fit = numpy.polyfit(fit_data_x,fit_data_y,1)
ccsdFit = exp(fit[1]), fit[0]

fit_data_x = [log(x) for x,y in beN_t_timing.data[-10:]]
fit_data_y = [log(y) for x,y in beN_t_timing.data[-10:]]
fit = numpy.polyfit(fit_data_x,fit_data_y,1)
_t_Fit = exp(fit[1]), fit[0]


def monomial(c,n):
	def f(x):  return  c * x**n
	return f


figure_0 = plot2d(ranges=((9.5,30),(0,360)), area=(300,200)).drawFrame()
figure_0.plotFunc(monomial(*ccsdFit),   connectors=blue)
figure_0.plotFunc(monomial(*_t_Fit), connectors=green)
figure_0.plotFunc(monomial(*XccsdFit),  connectors=tiger)
figure_0.plotData(beN_ccsd_timing.data,   connectors=False,markers=blue)
figure_0.plotData(beN_t_timing.data, connectors=False,markers=green)
figure_0.plotData(beN_Xccsd_timing.data,  connectors=False,markers=tiger)
figure_0.setXtics(periodicity=10, textform=text("{:.0f}",alignment=(center,top),  translate=(0,-10)))
figure_0.setYtics(periodicity=60,  textform=(text("{:.0f}",alignment=(right,center),translate=(-10,0)),to_min))
#figure_0 += text("Timings for Be$_n$", translate=(150,210),alignment=(center,bottom))
#figure_0 += text("min", translate=(-50,100),alignment=(right,center))
#figure_0 += text("$n$", translate=(150,-35),alignment=(center,top))
figure_0.pdf("timings",texlabels="texlabels/0/text*.svg")

print("E-CCSD:   {} * N**{}".format(*XccsdFit))
print("CCSD:     {} * N**{}".format(*ccsdFit))
print("(T):      {} * N**{}".format(*_t_Fit))


figure_4 = plot2d(ranges=((0,30),(-4e-4,0)), area=(300,200)).drawFrame()
figure_4.plotData(beN_ccsd_welldepth.data,  connectors=blue,markers=blue)
figure_4.plotData(beN_ccsd_t_welldepth.data,  connectors=green,markers=green)
figure_4.plotData(beN_Xccsd_welldepth.data,  connectors=tiger,markers=tiger)
figure_4.setXtics(periodicity=10, textform=text("{:.0f}",alignment=(center,top),  translate=(0,-10)))
figure_4.setYtics(periodicity=1e-4,  textform=(text("{:.2f}",alignment=(right,center),translate=(-10,0)),to_mEh))
figure_4.pdf("welldepth",texlabels="texlabels/4/text*.svg")
print("E-CCSD: {}%".format(100*beN_Xccsd_welldepth.data[28][1]/beN_ccsd_t_welldepth.data[29][1]))
print("CCSD:   {}%".format(100*beN_ccsd_welldepth.data[29][1] /beN_ccsd_t_welldepth.data[29][1]))
