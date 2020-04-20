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
Be_01.out:      * CCSD(T) total energy                  =  -14.612738068009943
Be_02.out:      * CCSD(T) total energy                  =  -29.225764276011677
Be_03.out:      * CCSD(T) total energy                  =  -43.838803002378697
Be_04.out:      * CCSD(T) total energy                  =  -58.451842945686359
Be_05.out:      * CCSD(T) total energy                  =  -73.064883066650822
Be_06.out:      * CCSD(T) total energy                  =  -87.677923230288684
Be_07.out:      * CCSD(T) total energy                  = -102.290963406314901
Be_08.out:      * CCSD(T) total energy                  = -116.904003586901027
Be_09.out:      * CCSD(T) total energy                  = -131.517043769511361
Be_10.out:      * CCSD(T) total energy                  = -146.130083953078554
Be_11.out:      * CCSD(T) total energy                  = -160.743124137100779
Be_12.out:      * CCSD(T) total energy                  = -175.356164321343385
Be_13.out:      * CCSD(T) total energy                  = -189.969204505695757
Be_14.out:      * CCSD(T) total energy                  = -204.582244690105767
Be_15.out:      * CCSD(T) total energy                  = -219.195284874547099
Be_16.out:      * CCSD(T) total energy                  = -233.808325059004119
Be_17.out:      * CCSD(T) total energy                  = -248.421365243469097
Be_18.out:      * CCSD(T) total energy                  = -263.034405427937770
Be_19.out:      * CCSD(T) total energy                  = -277.647445612406557
Be_20.out:      * CCSD(T) total energy                  = -292.260485796873752
Be_21.out:      * CCSD(T) total energy                  = -306.873525981337707
Be_22.out:      * CCSD(T) total energy                  = -321.486566165798934
Be_23.out:      * CCSD(T) total energy                  = -336.099606350254817
Be_24.out:      * CCSD(T) total energy                  = -350.712646534705868
Be_25.out:      * CCSD(T) total energy                  = -365.325686719153282
Be_26.out:      * CCSD(T) total energy                  = -379.938726903595295
Be_27.out:      * CCSD(T) total energy                  = -394.551767088031454
Be_28.out:      * CCSD(T) total energy                  = -409.164807272463463
Be_29.out:      * CCSD(T) total energy                  = -423.777847456890470
Be_30.out:      * CCSD(T) total energy                  = -438.390887641312077
""")
