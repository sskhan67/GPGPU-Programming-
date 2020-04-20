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
	dist   = float(fields[0].split("_")[2])
	energy = float(fields[-1])
	return dist, energy-2*atom.energy

data = X(
parse,
"""\
ccsd_Be2_3.0_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2174484874
ccsd_Be2_3.1_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2189912305
ccsd_Be2_3.2_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2202657983
ccsd_Be2_3.3_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2213163564
ccsd_Be2_3.4_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2221797804
ccsd_Be2_3.5_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2228860466
ccsd_Be2_3.6_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2234606339
ccsd_Be2_3.7_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2239252368
ccsd_Be2_3.8_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2242985814
ccsd_Be2_3.9_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2245965275
ccsd_Be2_4.0_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.224832393
ccsd_Be2_4.1_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2250174588
ccsd_Be2_4.2_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2251612455
ccsd_Be2_4.3_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2252717275
ccsd_Be2_4.4_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2253555359
ccsd_Be2_4.5_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254184801
ccsd_Be2_4.6_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254643326
ccsd_Be2_4.7_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.225497138
ccsd_Be2_4.8_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255198511
ccsd_Be2_4.9_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255348425
ccsd_Be2_5.0_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255439995
ccsd_Be2_5.1_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.225548813
ccsd_Be2_5.2_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255504503
ccsd_Be2_5.3_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255498174
ccsd_Be2_5.4_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255475524
ccsd_Be2_5.5_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255443129
ccsd_Be2_5.6_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255404228
ccsd_Be2_5.7_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255361758
ccsd_Be2_5.8_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255317851
ccsd_Be2_5.9_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255274023
ccsd_Be2_6.0_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255231316
ccsd_Be2_6.1_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.225519042
ccsd_Be2_6.2_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255151758
ccsd_Be2_6.3_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255115562
ccsd_Be2_6.4_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255081926
ccsd_Be2_6.5_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255050844
ccsd_Be2_6.6_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2255022246
ccsd_Be2_6.7_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254996019
ccsd_Be2_6.8_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254972023
ccsd_Be2_6.9_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254950106
ccsd_Be2_7.0_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254930111
ccsd_Be2_7.1_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254911882
ccsd_Be2_7.2_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254895267
ccsd_Be2_7.3_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254880125
ccsd_Be2_7.4_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254866322
ccsd_Be2_7.5_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254853734
ccsd_Be2_7.6_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254842248
ccsd_Be2_7.7_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254831761
ccsd_Be2_7.8_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254822177
ccsd_Be2_7.9_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254813412
ccsd_Be2_8.0_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254805389
ccsd_Be2_8.1_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254798038
ccsd_Be2_8.2_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254791297
ccsd_Be2_8.3_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254785108
ccsd_Be2_8.4_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254779421
ccsd_Be2_8.5_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254774191
ccsd_Be2_8.6_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254769375
ccsd_Be2_8.7_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254764938
ccsd_Be2_8.8_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254760844
ccsd_Be2_8.9_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254757065
ccsd_Be2_9.0_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254753573
ccsd_Be2_9.1_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254750342
ccsd_Be2_9.2_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254747352
ccsd_Be2_9.3_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.225474458
ccsd_Be2_9.4_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.225474201
ccsd_Be2_9.5_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254739625
ccsd_Be2_9.6_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254737409
ccsd_Be2_9.7_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254735348
ccsd_Be2_9.8_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254733431
ccsd_Be2_9.9_path=..:..:..:..:dimer_H:6-31g:H:0-115-0:load=H:0-115-0:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2254731645
""")
