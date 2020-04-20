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
ccsd_Be2_3.0_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2151409134
ccsd_Be2_3.1_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2171243195
ccsd_Be2_3.2_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2187464391
ccsd_Be2_3.3_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.220070924
ccsd_Be2_3.4_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2211506736
ccsd_Be2_3.5_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2220286384
ccsd_Be2_3.6_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2227404786
ccsd_Be2_3.7_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2233157519
ccsd_Be2_3.8_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2237791738
ccsd_Be2_3.9_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2241511335
ccsd_Be2_4.0_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2244483581
ccsd_Be2_4.1_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.224684695
ccsd_Be2_4.2_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2248716086
ccsd_Be2_4.3_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2250185632
ccsd_Be2_4.4_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2251333494
ccsd_Be2_4.5_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2252226999
ccsd_Be2_4.6_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2252911202
ccsd_Be2_4.7_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2253432704
ccsd_Be2_4.8_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2253826138
ccsd_Be2_4.9_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254119495
ccsd_Be2_5.0_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254335279
ccsd_Be2_5.1_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254491478
ccsd_Be2_5.2_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254602374
ccsd_Be2_5.3_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254679226
ccsd_Be2_5.4_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254730255
ccsd_Be2_5.5_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254763573
ccsd_Be2_5.6_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254783675
ccsd_Be2_5.7_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254794525
ccsd_Be2_5.8_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254799089
ccsd_Be2_5.9_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254799548
ccsd_Be2_6.0_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254797482
ccsd_Be2_6.1_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254794003
ccsd_Be2_6.2_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254789878
ccsd_Be2_6.3_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254785611
ccsd_Be2_6.4_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254781521
ccsd_Be2_6.5_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.225477779
ccsd_Be2_6.6_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254774505
ccsd_Be2_6.7_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254771693
ccsd_Be2_6.8_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254769341
ccsd_Be2_6.9_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254767412
ccsd_Be2_7.0_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254765857
ccsd_Be2_7.1_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254764624
ccsd_Be2_7.2_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.225476366
ccsd_Be2_7.3_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254762918
ccsd_Be2_7.4_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254762355
ccsd_Be2_7.5_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761934
ccsd_Be2_7.6_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761624
ccsd_Be2_7.7_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.22547614
ccsd_Be2_7.8_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761241
ccsd_Be2_7.9_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761132
ccsd_Be2_8.0_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761059
ccsd_Be2_8.1_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761012
ccsd_Be2_8.2_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254760985
ccsd_Be2_8.3_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254760971
ccsd_Be2_8.4_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254760967
ccsd_Be2_8.5_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254760969
ccsd_Be2_8.6_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254760976
ccsd_Be2_8.7_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254760984
ccsd_Be2_8.8_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254760995
ccsd_Be2_8.9_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761006
ccsd_Be2_9.0_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761016
ccsd_Be2_9.1_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761027
ccsd_Be2_9.2_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761037
ccsd_Be2_9.3_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761046
ccsd_Be2_9.4_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761054
ccsd_Be2_9.5_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761061
ccsd_Be2_9.6_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761068
ccsd_Be2_9.7_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761073
ccsd_Be2_9.8_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.2254761077
ccsd_Be2_9.9_path=:morescratch:adutoi:programs:Qode:Applications:Be_n:dimer_H:6-31g:run:molecule_CI:16-1-1:load=u_4.5_16-1-1_t_1e-6.pickle.log:Total Excitonic CCSD Energy =  -29.225476108
""")
