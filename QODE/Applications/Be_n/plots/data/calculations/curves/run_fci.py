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
from qode.util import output, textlog
import hamiltonian
import excitonic



# Usage:  python3 run_fci.py <n_atoms> <separation> <path>
n_atoms     = sys.argv[1]
separation  = sys.argv[2]
path        = sys.argv[3]

out = output(log=textlog(echo=True))

dist_mat = hamiltonian.dist_mat(float(separation), int(n_atoms))
H = hamiltonian.load_subH(path+"/{0}/{1}", dist_mat)

E_n = hamiltonian.E_nuc(dist_mat, 4)
E_e = excitonic.fci(H, out)

out.log("\nTotal Excitonic FCI Energy = ", E_n + E_e)
out.log.write("fci_Be{}_{}_path={}.log".format(n_atoms, separation, path.replace("/",":")))
