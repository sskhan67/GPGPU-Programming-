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
from data.extract import X
from data import atom

def parse(n, fields):
	dist   = float(".".join(fields[0].split("_")[1].split(".")[:2]))
	energy = float(fields[-1])
	return dist, energy-4*atom.energy

data = X(
parse,
"""\
Be4_3.8.out:      * CCSD(T) total energy                  =  -58.450155547197660
Be4_3.9.out:      * CCSD(T) total energy                  =  -58.450645128468203
Be4_4.0.out:      * CCSD(T) total energy                  =  -58.451029055338950
Be4_4.1.out:      * CCSD(T) total energy                  =  -58.451321862707132
Be4_4.2.out:      * CCSD(T) total energy                  =  -58.451537468728588
Be4_4.3.out:      * CCSD(T) total energy                  =  -58.451688718831242
Be4_4.4.out:      * CCSD(T) total energy                  =  -58.451787155632850
Be4_4.5.out:      * CCSD(T) total energy                  =  -58.451842942707174
Be4_4.6.out:      * CCSD(T) total energy                  =  -58.451864887450000
Be4_4.7.out:      * CCSD(T) total energy                  =  -58.451860523913659
Be4_4.8.out:      * CCSD(T) total energy                  =  -58.451836226809533
Be4_4.9.out:      * CCSD(T) total energy                  =  -58.451797337066097
Be4_5.0.out:      * CCSD(T) total energy                  =  -58.451748286239891
Be4_5.1.out:      * CCSD(T) total energy                  =  -58.451692712586308
Be4_5.2.out:      * CCSD(T) total energy                  =  -58.451633565427869
Be4_5.3.out:      * CCSD(T) total energy                  =  -58.451573197360695
Be4_5.4.out:      * CCSD(T) total energy                  =  -58.451513445157921
Be4_5.5.out:      * CCSD(T) total energy                  =  -58.451455701324875
Be4_5.6.out:      * CCSD(T) total energy                  =  -58.451400978000294
Be4_5.7.out:      * CCSD(T) total energy                  =  -58.451349964851460
Be4_5.8.out:      * CCSD(T) total energy                  =  -58.451303082442010
Be4_5.9.out:      * CCSD(T) total energy                  =  -58.451260529160798
Be4_6.0.out:      * CCSD(T) total energy                  =  -58.451222329246669
Be4_6.1.out:      * CCSD(T) total energy                  =  -58.451188369827818
Be4_6.2.out:      * CCSD(T) total energy                  =  -58.451158438437481
Be4_6.3.out:      * CCSD(T) total energy                  =  -58.451132254096947
Be4_6.4.out:      * CCSD(T) total energy                  =  -58.451109493800466
Be4_6.5.out:      * CCSD(T) total energy                  =  -58.451089814484710
Be4_6.6.out:      * CCSD(T) total energy                  =  -58.451072869924587
Be4_6.7.out:      * CCSD(T) total energy                  =  -58.451058324455282
Be4_6.8.out:      * CCSD(T) total energy                  =  -58.451045861622468
Be4_6.9.out:      * CCSD(T) total energy                  =  -58.451035190409662
Be4_7.0.out:      * CCSD(T) total energy                  =  -58.451026048874098
Be4_7.1.out:      * CCSD(T) total energy                  =  -58.451018205067193
Be4_7.2.out:      * CCSD(T) total energy                  =  -58.451011456748489
Be4_7.3.out:      * CCSD(T) total energy                  =  -58.451005630035120
Be4_7.4.out:      * CCSD(T) total energy                  =  -58.451000576993785
Be4_7.5.out:      * CCSD(T) total energy                  =  -58.450996173036771
Be4_7.6.out:      * CCSD(T) total energy                  =  -58.450992314116384
Be4_7.7.out:      * CCSD(T) total energy                  =  -58.450988913886242
Be4_7.8.out:      * CCSD(T) total energy                  =  -58.450985901070368
Be4_7.9.out:      * CCSD(T) total energy                  =  -58.450983216990849
Be4_8.0.out:      * CCSD(T) total energy                  =  -58.450980813463332
Be4_8.1.out:      * CCSD(T) total energy                  =  -58.450978650878668
Be4_8.2.out:      * CCSD(T) total energy                  =  -58.450976696664348
Be4_8.3.out:      * CCSD(T) total energy                  =  -58.450974923851284
Be4_8.4.out:      * CCSD(T) total energy                  =  -58.450973310082183
Be4_8.5.out:      * CCSD(T) total energy                  =  -58.450971836662568
Be4_8.6.out:      * CCSD(T) total energy                  =  -58.450970487850718
Be4_8.7.out:      * CCSD(T) total energy                  =  -58.450969250290591
Be4_8.8.out:      * CCSD(T) total energy                  =  -58.450968112541439
Be4_8.9.out:      * CCSD(T) total energy                  =  -58.450967064731238
Be4_9.0.out:      * CCSD(T) total energy                  =  -58.450966098269760
Be4_9.1.out:      * CCSD(T) total energy                  =  -58.450965205624115
Be4_9.2.out:      * CCSD(T) total energy                  =  -58.450964380144413
Be4_9.3.out:      * CCSD(T) total energy                  =  -58.450963615925318
Be4_9.4.out:      * CCSD(T) total energy                  =  -58.450962907694112
Be4_9.5.out:      * CCSD(T) total energy                  =  -58.450962250722448
Be4_9.6.out:      * CCSD(T) total energy                  =  -58.450961640751586
Be4_9.7.out:      * CCSD(T) total energy                  =  -58.450961073937876
Be4_9.8.out:      * CCSD(T) total energy                  =  -58.450960546799259
Be4_9.9.out:      * CCSD(T) total energy                  =  -58.450960056174488
""")
