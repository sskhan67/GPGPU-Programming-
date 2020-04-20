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
# utilities
import sys
import time
from qode.util import parallel, output, textlog, indent

# SCF and CCSD modules
from qode.SCF_diag_full                 import hf as scf_routine
import ccsd
import numpy as np


np.set_printoptions(threshold=np.nan,linewidth=270, precision=10)

out = output(log=textlog(echo=True))
resources = parallel.resources(int(sys.argv[1]))
input_file = 'h2'
t1 = time.time()

out.log('\n### HARTREE-FOCK CALCULATION ###\n')
E_0, KE, NAE, NRE, ERE, F_mat, T_mat, N_mat, V_mat, vHF, num_alpha_elec, num_beta_elec, num_spin_orb = scf_routine.main(input_file+'.in', 1e-6)


np.save( input_file + "_h_mat.npy", T_mat + N_mat )
np.save( input_file + "_V_mat.npy", V_mat )

# print(V_mat)

# out.log("Parsing CC input ...")
h_mat = T_mat + N_mat
print("h mat:")
print(h_mat)

# occ_orbs = []
# vrt_orbs = []
# for p in range(num_alpha_elec):  occ_orbs += [p]
# for p in range(num_beta_elec):   occ_orbs += [p+(num_spin_orb//2)]
# for p in range(num_spin_orb):
# 	if p not in occ_orbs:  vrt_orbs += [p]
# out.log(indent("occupied orbital indices = ", occ_orbs))
# out.log(indent("virtual  orbital indices = ", vrt_orbs))

#out.log('### COUPLED-CLUSTER CALCULATION ###\n')
#out.ccsd = ccsd.main(h_mat, V_mat, F_mat, NRE, occ_orbs, vrt_orbs, textlog=out.log.sublog(), resources=resources)

#t2 = time.time()
#out.log('\n### RESULTS ###\n')
#out.log("TOTAL ENERGY =", out.ccsd.energy)
#out.log("TOTAL TIME   =", t2 - t1, "seconds\n")
#out.log.write(input_file + '.out')
