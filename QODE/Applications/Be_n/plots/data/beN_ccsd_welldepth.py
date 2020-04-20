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
	num    = int(fields[0].split(".")[0].split("_")[1])
	energy = float(fields[-1])/num - atom.energy
	return num, energy

data = X(
parse,
"""\
Be_01.out:      * CCSD total energy                     =  -14.612738068009943
Be_02.out:      * CCSD total energy                     =  -29.225703143802811
Be_03.out:      * CCSD total energy                     =  -43.838675973455828
Be_04.out:      * CCSD total energy                     =  -58.451649645724331
Be_05.out:      * CCSD total energy                     =  -73.064623452540772
Be_06.out:      * CCSD total energy                     =  -87.677597294355223
Be_07.out:      * CCSD total energy                     = -102.290571146610390
Be_08.out:      * CCSD total energy                     = -116.903545002785592
Be_09.out:      * CCSD total energy                     = -131.516518860734067
Be_10.out:      * CCSD total energy                     = -146.129492719526496
Be_11.out:      * CCSD total energy                     = -160.742466578717313
Be_12.out:      * CCSD total energy                     = -175.355440438097588
Be_13.out:      * CCSD total energy                     = -189.968414297569694
Be_14.out:      * CCSD total energy                     = -204.581388157088554
Be_15.out:      * CCSD total energy                     = -219.194362016631828
Be_16.out:      * CCSD total energy                     = -233.807335876186329
Be_17.out:      * CCSD total energy                     = -248.420309735745718
Be_18.out:      * CCSD total energy                     = -263.033283595306671
Be_19.out:      * CCSD total energy                     = -277.646257454866372
Be_20.out:      * CCSD total energy                     = -292.259231314423346
Be_21.out:      * CCSD total energy                     = -306.872205173976340
Be_22.out:      * CCSD total energy                     = -321.485179033526038
Be_23.out:      * CCSD total energy                     = -336.098152893069937
Be_24.out:      * CCSD total energy                     = -350.711126752608664
Be_25.out:      * CCSD total energy                     = -365.324100612143468
Be_26.out:      * CCSD total energy                     = -379.937074471672759
Be_27.out:      * CCSD total energy                     = -394.550048331195967
Be_28.out:      * CCSD total energy                     = -409.163022190714912
Be_29.out:      * CCSD total energy                     = -423.775996050228741
Be_30.out:      * CCSD total energy                     = -438.388969909737114
""")
