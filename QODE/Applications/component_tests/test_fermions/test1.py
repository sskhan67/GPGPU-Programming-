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
from qode.fermion_field import config, occ_strings



occ_orbs = [0, 4]
vrt_orbs = [1,2,3, 5,6,7]




print("FCI occupation strings for 2 electrons in 4 spin orbitals")
strings = occ_strings.all_occ_strings(4,2)
for string in strings:  print(strings.index(string),string)
print("LENGTH =",len(strings))

print("CISD occupation strings for 2 electrons in 8 spin orbitals blocked by alpha and beta")
cisd_strings = occ_strings.CISD(occ_orbs,vrt_orbs)
for string in cisd_strings:  print(cisd_strings.index(string),string)
print("LENGTH =",len(cisd_strings))




"""\
Here we really just test the __init__ function and the print function by conjuring up a bunch of occupancy strings
and wrapping them up in configuration objects.

The configurations mimuc CISD on a 2e- system (so FCI) where there are 4 spatial orbitals blocked into alpha and beta spin.
"""
cisd_configs = occ_strings.CISD(occ_orbs,vrt_orbs)
for configuration in cisd_configs:  print(config.configuration(configuration))
print("LENGTH =",len(cisd_configs))
