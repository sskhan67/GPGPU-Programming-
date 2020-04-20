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
from qode.fermion_field import occ_strings
from qode.fermion_field.state import *

all_occ_strings = occ_strings.all_occ_strings(8,4)

print(state(all_occ_strings,0))

n_states = len(all_occ_strings)
cum = 0.
for i in range(n_states):
	for j in range(n_states):
		inner_pdt = dot(state(all_occ_strings,i),state(all_occ_strings,j))
		if i==j:  cum += (inner_pdt - 1)**2
		else:     cum +=  inner_pdt**2
print(cum,"\n")

occ_orbs = [0, 4]
vrt_orbs = [1,2,3, 5,6,7]
cisd_strings = occ_strings.CISD(occ_orbs,vrt_orbs)
print(state(cisd_strings,0))
