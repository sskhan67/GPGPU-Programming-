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
import time
from qode.util import output, indent

# operators, SCF and CCSD modules
from qode                               import fermion_field
from qode.many_body.nested_operator     import operator, space_traits, BakerCampbellHausdorff
from qode.many_body                     import CCSD

# Hamiltonian
from hamiltonian import hamiltonian, inv_fock



def main(h_mat, V_mat, F_mat, E_nuclear_repul, occ_orbs, vrt_orbs, textlog, resources):
	out = output(log=textlog)

	out.log("Partioning fluctuations (excitations, flat, deexcitations) ...")
	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)

	out.log("Building normal-ordered H ...")
	H = hamiltonian(E_nuclear_repul, h_mat, V_mat, fluctuations)

	out.log("Getting orbital energy preconditioner ...")
	invF = inv_fock(F_mat, fluctuations)

	out.log("Initializing T ...")
	T = operator(fluctuations)		# Initialized to zero

	BCH = BakerCampbellHausdorff(fluctuations, resources)
	out.log("Baker-Campbell-Hausdorff module initialized.")

	out.log("Perform coupled-cluster iterations.")
	t_ccsd_start = time.time()
	E = CCSD(H, T, invF, BCH, space_traits, out.log.sublog())
	t_ccsd_end = time.time()

	out.log("CCSD Energy =", E, 'Hartree')
	out.log("CCSD TIME {}".format(t_ccsd_end - t_ccsd_start))

	return out(energy=E)
