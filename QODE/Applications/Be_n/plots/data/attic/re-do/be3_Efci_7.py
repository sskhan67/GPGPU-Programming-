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
fci_Be3_4.3_path=..:..:dimer_H:6-31g:run:molecule_CI:16-100-100:load=u_4.5_16-100-100_t_1e-7.pickle.log:Total Excitonic FCI Energy =  -43.8386989625
fci_Be3_4.5_path=..:..:dimer_H:6-31g:run:molecule_CI:16-100-100:load=u_4.5_16-100-100_t_1e-7.pickle.log:Total Excitonic FCI Energy =  -43.8388247002
fci_Be3_4.7_path=..:..:dimer_H:6-31g:run:molecule_CI:16-100-100:load=u_4.5_16-100-100_t_1e-7.pickle.log:Total Excitonic FCI Energy =  -43.8388430094
""")
