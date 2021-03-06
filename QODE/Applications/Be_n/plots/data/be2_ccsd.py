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
	return dist, energy-2*atom.energy

data = X(
parse,
"""\
Be2_3.8.out:      * CCSD total energy                     =  -29.224967046654065
Be2_3.9.out:      * CCSD total energy                     =  -29.225171186615817
Be2_4.0.out:      * CCSD total energy                     =  -29.225332903176852
Be2_4.1.out:      * CCSD total energy                     =  -29.225458346733205
Be2_4.2.out:      * CCSD total energy                     =  -29.225553228132895
Be2_4.3.out:      * CCSD total energy                     =  -29.225622708973948
Be2_4.4.out:      * CCSD total energy                     =  -29.225671356353036
Be2_4.5.out:      * CCSD total energy                     =  -29.225703143181086
Be2_4.6.out:      * CCSD total energy                     =  -29.225721478435492
Be2_4.7.out:      * CCSD total energy                     =  -29.225729254616130
Be2_4.8.out:      * CCSD total energy                     =  -29.225728903482548
Be2_4.9.out:      * CCSD total energy                     =  -29.225722452648753
Be2_5.0.out:      * CCSD total energy                     =  -29.225711580071700
Be2_5.1.out:      * CCSD total energy                     =  -29.225697662653264
Be2_5.2.out:      * CCSD total energy                     =  -29.225681819426839
Be2_5.3.out:      * CCSD total energy                     =  -29.225664948373236
Be2_5.4.out:      * CCSD total energy                     =  -29.225647757977047
Be2_5.5.out:      * CCSD total energy                     =  -29.225630794478544
Be2_5.6.out:      * CCSD total energy                     =  -29.225614465644671
Be2_5.7.out:      * CCSD total energy                     =  -29.225599060973128
Be2_5.8.out:      * CCSD total energy                     =  -29.225584772561401
Be2_5.9.out:      * CCSD total energy                     =  -29.225571710058130
Be2_6.0.out:      * CCSD total energy                     =  -29.225559918564002
Be2_6.1.out:      * CCSD total energy                     =  -29.225549391199038
Be2_6.2.out:      * CCSD total energy                     =  -29.225540083273071
Be2_6.3.out:      * CCSD total energy                     =  -29.225531922360414
Be2_6.4.out:      * CCSD total energy                     =  -29.225524818522231
Be2_6.5.out:      * CCSD total energy                     =  -29.225518671853056
Be2_6.6.out:      * CCSD total energy                     =  -29.225513378927353
Be2_6.7.out:      * CCSD total energy                     =  -29.225508837452114
Be2_6.8.out:      * CCSD total energy                     =  -29.225504949873027
Be2_6.9.out:      * CCSD total energy                     =  -29.225501625697060
Be2_7.0.out:      * CCSD total energy                     =  -29.225498782818022
Be2_7.1.out:      * CCSD total energy                     =  -29.225496348314710
Be2_7.2.out:      * CCSD total energy                     =  -29.225494258317131
Be2_7.3.out:      * CCSD total energy                     =  -29.225492457791322
Be2_7.4.out:      * CCSD total energy                     =  -29.225490899900002
Be2_7.5.out:      * CCSD total energy                     =  -29.225489545108285
Be2_7.6.out:      * CCSD total energy                     =  -29.225488360467399
Be2_7.7.out:      * CCSD total energy                     =  -29.225487318610043
Be2_7.8.out:      * CCSD total energy                     =  -29.225486397001688
Be2_7.9.out:      * CCSD total energy                     =  -29.225485577126292
Be2_8.0.out:      * CCSD total energy                     =  -29.225484843807813
Be2_8.1.out:      * CCSD total energy                     =  -29.225484184615937
Be2_8.2.out:      * CCSD total energy                     =  -29.225483589339781
Be2_8.3.out:      * CCSD total energy                     =  -29.225483049584973
Be2_8.4.out:      * CCSD total energy                     =  -29.225482558402778
Be2_8.5.out:      * CCSD total energy                     =  -29.225482110006951
Be2_8.6.out:      * CCSD total energy                     =  -29.225481699548205
Be2_8.7.out:      * CCSD total energy                     =  -29.225481322917968
Be2_8.8.out:      * CCSD total energy                     =  -29.225480976614733
Be2_8.9.out:      * CCSD total energy                     =  -29.225480657624978
Be2_9.0.out:      * CCSD total energy                     =  -29.225480363305422
Be2_9.1.out:      * CCSD total energy                     =  -29.225480091440662
Be2_9.2.out:      * CCSD total energy                     =  -29.225479839937670
Be2_9.3.out:      * CCSD total energy                     =  -29.225479607030479
Be2_9.4.out:      * CCSD total energy                     =  -29.225479391121358
Be2_9.5.out:      * CCSD total energy                     =  -29.225479190776348
Be2_9.6.out:      * CCSD total energy                     =  -29.225479004711552
Be2_9.7.out:      * CCSD total energy                     =  -29.225478831759816
Be2_9.8.out:      * CCSD total energy                     =  -29.225478670867254
Be2_9.9.out:      * CCSD total energy                     =  -29.225478521076923
""")
