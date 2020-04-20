from pytoon        import *
from pytoon.macros import *

tiger = "#F47920"
tiger_cross = composite()
tiger_cross += line(begin=(-0.5,-0.5),end=(0.5,0.5),color=tiger,weight=1/3.)
tiger_cross += line(begin=(0.5,-0.5),end=(-0.5,0.5),color=tiger,weight=1/3.)
tiger_cross = tiger_cross(scale=3)

textX = text(xheight=15)

def to_mEh(nrg):  return 1000*nrg



from data import fci
from data import xr_20
from data import xr_30
from data import xr_35
from data import xr_40
from data import xr_45
from data import xr_100

dimer = plot2d(ranges=((2.99,7.99),(-5e-4,5e-4)), area=(300,200)).drawFrame()

dimer.plotData(fci.data,    connectors=(black,3), markers=False)
dimer.plotData(xr_20.data,  connectors=False,     markers=circle_pt(fill=orange,scale=5))
dimer.plotData(xr_30.data,  connectors=False,     markers=circle_pt(fill=magenta,scale=5))
dimer.plotData(xr_35.data,  connectors=False,     markers=circle_pt(fill=blue,scale=5))
dimer.plotData(xr_40.data,  connectors=False,     markers=circle_pt(fill=green,scale=5))
dimer.plotData(xr_45.data,  connectors=False,     markers=circle_pt(fill=red,scale=5))
dimer.plotData(xr_100.data, connectors=False,     markers=tiger_cross)

dimer.setXtics(periodicity=1,   textform=textX("{:.0f}",alignment=(center,top),  translate=(0,-10)))
dimer.setYtics(periodicity=2e-3,textform=(textX("{:.1f}",alignment=(right,center),translate=(-10,-0)),to_mEh))
dimer.pdf("PE",texlabels="texlabels/PE/text*.svg")
