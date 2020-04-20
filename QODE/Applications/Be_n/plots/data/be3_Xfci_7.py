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
	return dist, energy-3*atom.energy

data = X(
parse,
"""\
fci_Be3_3.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8282727889
fci_Be3_3.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8299585721
fci_Be3_3.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8315436114
fci_Be3_3.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8329721307
fci_Be3_3.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8342173854
fci_Be3_3.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8352724166
fci_Be3_3.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8361449372
fci_Be3_3.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.836851118
fci_Be3_3.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8374115575
fci_Be3_3.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8378476134
fci_Be3_4.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8381794667
fci_Be3_4.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.83842552
fci_Be3_4.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8386019875
fci_Be3_4.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8387227659
fci_Be3_4.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8387995391
fci_Be3_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838842668
fci_Be3_4.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8388586627
fci_Be3_4.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8388547852
fci_Be3_4.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838836282
fci_Be3_4.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838807379
fci_Be3_5.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8387714621
fci_Be3_5.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8387312257
fci_Be3_5.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8386887937
fci_Be3_5.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8386458178
fci_Be3_5.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8386034407
fci_Be3_5.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8385628541
fci_Be3_5.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8385245698
fci_Be3_5.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8384890221
fci_Be3_5.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8384564605
fci_Be3_5.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8384269857
fci_Be3_6.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838400582
fci_Be3_6.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383771456
fci_Be3_6.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383565094
fci_Be3_6.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383384646
fci_Be3_6.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383227777
fci_Be3_6.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383092055
fci_Be3_6.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382975058
fci_Be3_6.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382874457
fci_Be3_6.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382788076
fci_Be3_6.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382713924
fci_Be3_7.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382650215
fci_Be3_7.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382595374
fci_Be3_7.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382548032
fci_Be3_7.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382507011
fci_Be3_7.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382471312
fci_Be3_7.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382440091
fci_Be3_7.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382412645
fci_Be3_7.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382388388
fci_Be3_7.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382366837
fci_Be3_7.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838234759
fci_Be3_8.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382330319
fci_Be3_8.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382314751
fci_Be3_8.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382300663
fci_Be3_8.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382287867
fci_Be3_8.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382276206
fci_Be3_8.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382262599
fci_Be3_8.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382252864
fci_Be3_8.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838224393
fci_Be3_8.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382235715
fci_Be3_8.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382228149
fci_Be3_9.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382221169
fci_Be3_9.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382214723
fci_Be3_9.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382208761
fci_Be3_9.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382203242
fci_Be3_9.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382198127
fci_Be3_9.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382193382
fci_Be3_9.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382188976
fci_Be3_9.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382184881
fci_Be3_9.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382181073
fci_Be3_9.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-7:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382177527
""")
