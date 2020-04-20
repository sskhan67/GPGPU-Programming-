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
	return dist, energy-3*atom.energy

data = X(
parse,
"""\
Be3_3.8.out:      * CCSD total energy                     =  -43.837229084077052
Be3_3.9.out:      * CCSD total energy                     =  -43.837629466217464
Be3_4.0.out:      * CCSD total energy                     =  -43.837947149515969
Be3_4.1.out:      * CCSD total energy                     =  -43.838193922415059
Be3_4.2.out:      * CCSD total energy                     =  -43.838380775004133
Be3_4.3.out:      * CCSD total energy                     =  -43.838517679774583
Be3_4.4.out:      * CCSD total energy                     =  -43.838613497818606
Be3_4.5.out:      * CCSD total energy                     =  -43.838675972253171
Be3_4.6.out:      * CCSD total energy                     =  -43.838711778203049
Be3_4.7.out:      * CCSD total energy                     =  -43.838726606480684
Be3_4.8.out:      * CCSD total energy                     =  -43.838725263137150
Be3_4.9.out:      * CCSD total energy                     =  -43.838711772651301
Be3_5.0.out:      * CCSD total energy                     =  -43.838689477361044
Be3_5.1.out:      * CCSD total energy                     =  -43.838661128203967
Be3_5.2.out:      * CCSD total energy                     =  -43.838628965247466
Be3_5.3.out:      * CCSD total energy                     =  -43.838594787891559
Be3_5.4.out:      * CCSD total energy                     =  -43.838560015957711
Be3_5.5.out:      * CCSD total energy                     =  -43.838525742883640
Be3_5.6.out:      * CCSD total energy                     =  -43.838492783029523
Be3_5.7.out:      * CCSD total energy                     =  -43.838461713857903
Be3_5.8.out:      * CCSD total energy                     =  -43.838432914606486
Be3_5.9.out:      * CCSD total energy                     =  -43.838406601233935
Be3_6.0.out:      * CCSD total energy                     =  -43.838382858696590
Be3_6.1.out:      * CCSD total energy                     =  -43.838361669098461
Be3_6.2.out:      * CCSD total energy                     =  -43.838342938962690
Be3_6.3.out:      * CCSD total energy                     =  -43.838326520157871
Be3_6.4.out:      * CCSD total energy                     =  -43.838312229793985
Be3_6.5.out:      * CCSD total energy                     =  -43.838299865447276
Be3_6.6.out:      * CCSD total energy                     =  -43.838289218359755
Be3_6.7.out:      * CCSD total energy                     =  -43.838280082221992
Be3_6.8.out:      * CCSD total energy                     =  -43.838272260636963
Be3_6.9.out:      * CCSD total energy                     =  -43.838265571515500
Be3_7.0.out:      * CCSD total energy                     =  -43.838259849845620
Be3_7.1.out:      * CCSD total energy                     =  -43.838254949015230
Be3_7.2.out:      * CCSD total energy                     =  -43.838250740719218
Be3_7.3.out:      * CCSD total energy                     =  -43.838247114411693
Be3_7.4.out:      * CCSD total energy                     =  -43.838243975981953
Be3_7.5.out:      * CCSD total energy                     =  -43.838241246094661
Be3_7.6.out:      * CCSD total energy                     =  -43.838238858495657
Be3_7.7.out:      * CCSD total energy                     =  -43.838236758277056
Be3_7.8.out:      * CCSD total energy                     =  -43.838234900120860
Be3_7.9.out:      * CCSD total energy                     =  -43.838233246813274
Be3_8.0.out:      * CCSD total energy                     =  -43.838231767848804
Be3_8.1.out:      * CCSD total energy                     =  -43.838230438261242
Be3_8.2.out:      * CCSD total energy                     =  -43.838229237496989
Be3_8.3.out:      * CCSD total energy                     =  -43.838228148662068
Be3_8.4.out:      * CCSD total energy                     =  -43.838227157781482
Be3_8.5.out:      * CCSD total energy                     =  -43.838226253193312
Be3_8.6.out:      * CCSD total energy                     =  -43.838225425124399
Be3_8.7.out:      * CCSD total energy                     =  -43.838224665306477
Be3_8.8.out:      * CCSD total energy                     =  -43.838223966674413
Be3_8.9.out:      * CCSD total energy                     =  -43.838223323152548
Be3_9.0.out:      * CCSD total energy                     =  -43.838222729459716
Be3_9.1.out:      * CCSD total energy                     =  -43.838222180977191
Be3_9.2.out:      * CCSD total energy                     =  -43.838221673611443
Be3_9.3.out:      * CCSD total energy                     =  -43.838221203807755
Be3_9.4.out:      * CCSD total energy                     =  -43.838220768276173
Be3_9.5.out:      * CCSD total energy                     =  -43.838220364163959
Be3_9.6.out:      * CCSD total energy                     =  -43.838219988857780
Be3_9.7.out:      * CCSD total energy                     =  -43.838219640009456
Be3_9.8.out:      * CCSD total energy                     =  -43.838219315491813
Be3_9.9.out:      * CCSD total energy                     =  -43.838219013374157
""")
