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
fci_Be3_3.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8278940935
fci_Be3_3.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8297172446
fci_Be3_3.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8313902005
fci_Be3_3.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8328740639
fci_Be3_3.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8341540145
fci_Be3_3.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8352309486
fci_Be3_3.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8361175024
fci_Be3_3.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.836832876
fci_Be3_3.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8373994312
fci_Be3_3.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8378396101
fci_Be3_4.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838174248
fci_Be3_4.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8384221649
fci_Be3_4.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8385998439
fci_Be3_4.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8387213713
fci_Be3_4.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8387985666
fci_Be3_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8388418919
fci_Be3_4.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8388579326
fci_Be3_4.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8388540089
fci_Be3_4.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8388354089
fci_Be3_4.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8388063889
fci_Be3_5.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8387703553
fci_Be3_5.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8387300153
fci_Be3_5.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8386874996
fci_Be3_5.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8386444626
fci_Be3_5.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8386020465
fci_Be3_5.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8385614409
fci_Be3_5.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8385231545
fci_Be3_5.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8384876181
fci_Be3_5.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8384550775
fci_Be3_5.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8384256304
fci_Be3_6.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383992583
fci_Be3_6.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383758555
fci_Be3_6.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383552533
fci_Be3_6.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383372415
fci_Be3_6.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383215862
fci_Be3_6.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8383080434
fci_Be3_6.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382963708
fci_Be3_6.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382863355
fci_Be3_6.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382777198
fci_Be3_6.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382703248
fci_Be3_7.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838263972
fci_Be3_7.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382585042
fci_Be3_7.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382537845
fci_Be3_7.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382496954
fci_Be3_7.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382461372
fci_Be3_7.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382430256
fci_Be3_7.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382402904
fci_Be3_7.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382378732
fci_Be3_7.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382357257
fci_Be3_7.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838233808
fci_Be3_8.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382320872
fci_Be3_8.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382305362
fci_Be3_8.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382291325
fci_Be3_8.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382278577
fci_Be3_8.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838226696
fci_Be3_8.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382253391
fci_Be3_8.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382243693
fci_Be3_8.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382234793
fci_Be3_8.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382226609
fci_Be3_8.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382219071
fci_Be3_9.0_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382212118
fci_Be3_9.1_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382205695
fci_Be3_9.2_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382199756
fci_Be3_9.3_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382194258
fci_Be3_9.4_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382189162
fci_Be3_9.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382184435
fci_Be3_9.6_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382180045
fci_Be3_9.7_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382175966
fci_Be3_9.8_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.8382172172
fci_Be3_9.9_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic FCI Energy =  -43.838216864
""")
