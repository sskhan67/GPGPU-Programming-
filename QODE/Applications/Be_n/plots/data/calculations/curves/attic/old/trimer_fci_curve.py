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
import sys
from   qode.util import parallel, output, textlog
import hamiltonian
import fci

# Usage:  python3 trimer_fci_curve.py <basis> <state-num-string> <compression_tag>
basis       = sys.argv[1]
state_nums  = sys.argv[2]
compression = sys.argv[3]

filestem_pattern = "dimer_H/{basis}/run/molecule_CI/{{1}}_{{0}}_{state_nums}_{compression}".format(basis=basis,state_nums=state_nums,compression=compression)

out  = output(log=textlog(echo=False))

distances = [\
"4.0", "4.1", "4.2", "4.3", "4.4", "4.5", "4.6", "4.7", "4.8", "4.9", \
"5.0", "5.1", "5.2", "5.3", "5.4", "5.5", "5.6", "5.7", "5.8", "5.9", \
"6.0" ]

for dist in distances:
	doubled = "{:.1f}".format(2*float(dist))
	dist_mat = [[".",dist,doubled],[".",".",dist],[".",".","."]]
	E_e = fci.main(hamiltonian.load_subH(filestem_pattern, dist_mat), out)
	E   = E_e + hamiltonian.E_nuc(dist_mat,4)
	print("DATA:  ",dist,E)	

out.log.write("trimer_fci_curve.log")
