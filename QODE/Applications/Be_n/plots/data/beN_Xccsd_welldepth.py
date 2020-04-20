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
	num  = int(fields[0].split("_")[1][2:])
	energy = float(fields[-1])/num - atom.energy
	return num, energy

data = X(
parse,
"""\
ccsd_Be2_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -29.2258034218
ccsd_Be3_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -43.8387870429
ccsd_Be4_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -58.4517695208
ccsd_Be5_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -73.064751709
ccsd_Be6_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -87.6777338214
ccsd_Be7_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -102.290716104
ccsd_Be8_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -116.903698548
ccsd_Be9_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -131.516681151
ccsd_Be10_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -146.129663908
ccsd_Be11_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -160.742646818
ccsd_Be12_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -175.355629879
ccsd_Be13_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -189.968613091
ccsd_Be14_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -204.58159645
ccsd_Be15_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -219.194579958
ccsd_Be16_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -233.807563613
ccsd_Be17_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -248.420547414
ccsd_Be18_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -263.033531362
ccsd_Be19_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -277.646515455
ccsd_Be20_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -292.259499693
ccsd_Be21_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -306.872484076
ccsd_Be22_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -321.485468603
ccsd_Be23_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -336.098453274
ccsd_Be24_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -350.711438089
ccsd_Be25_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -365.324423048
ccsd_Be26_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -379.937408149
ccsd_Be27_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -394.550393394
ccsd_Be28_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -409.163378782
ccsd_Be29_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -423.776364313
ccsd_Be30_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -438.389349986
ccsd_Be40_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -584.519214527
ccsd_Be50_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -730.649093195
ccsd_Be60_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -876.778985918
ccsd_Be70_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -1022.90889265
ccsd_Be80_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -1169.03881336
ccsd_Be90_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -1315.16874803
ccsd_Be100_4.5_path=..:..:..:..:dimer_H:6-31g:H:16-115-550:load=H:16-115-550:thresh=1e-6:4.5:u.pickle.log:Total Excitonic CCSD Energy =  -1461.29869662
""")
