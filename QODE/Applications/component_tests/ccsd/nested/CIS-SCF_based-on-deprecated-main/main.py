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

from qode.math import linear_inner_product_space, lanczos
from qode.many_body.nested_operator import mask
from util import rotate_electronic_scf

# operators, SCF and CCSD modules
from qode                           import fermion_field
from qode.SCF_diag_full             import hfHacked as scf_routine_hacked
from qode.many_body.nested_operator import operator, space_traits, BakerCampbellHausdorff
from qode.many_body                 import CCSD

# Hamiltonian
from hamiltonian import hamiltonian, digest_V_mat, PT_step_conditioner, excite_mask, const_excite_projection, super_operator_wrapper		# move some of these to nested_operator module



def main():												# main function called at the end of this file
	out = output(log=textlog(echo=True))
	resources = parallel.resources(int(sys.argv[1]))

	t1 = time.time()

	integrals, orbitals = parse_input()
	textlog("Input parsed.")
	textlog(indent("occupied orbital indices = ", orbitals.occ))
	textlog(indent("virtual  orbital indices = ", orbitals.vrt))

	fluctuations = fermion_field.OV_partitioning(orbitals.occ, orbitals.vrt)
	textlog("Fluctuations partitioned (excitations, flat, deexcitations).")

	out.log("Performing generalized SCF algorithm ...")
	fluctuations, normal_ordered_H, SCF_energy_diffs = general_scf(fluctuations, integrals, rotate_electronic_scf(), textlog=out.log.sublog(), resources=resources)
	raise Exception("please be right?!")

	pseudoHessian = PT_step_conditioner(out.SCF.energy_diffs, fluctuations)
	out.log("Constant T-space pseudo-Hessian constructed from perturbation theory.")

	out.log("Performing generalized CC algorithm ...")
	out.CCSD = ccsd(fluctuations, normal_ordered_H, pseudoHessian, textlog=out.log.sublog(), resources=resources)

	t2 = time.time()

	out.log("Summary:")
	out.log("TOTAL ENERGY =", out.CCSD.energy)
	out.log("TOTAL TIME   =", t2 - t1, "seconds\n")
	out.log.write('1.0.out')



### Major functions

def general_scf(fluctuations, integrals, rotate_scf, textlog, resources, thresh=1e-4, krylov_dim=5):
	textlog("==================================================================================")
	Enuc, h1, V2 = integrals.Enuc, integrals.h1, integrals.V2
	textlog("Integrals unpacked.")
	#
	Cspace = linear_inner_product_space(space_traits)
	Pop  = Cspace.lin_op(const_excite_projection)
	textlog("CIRS space set up (CIS, including reference).")
	#
	iteration = 0
	rot_norm = float("inf")
	while rot_norm>thresh:
		textlog("==================================================================================")
		textlog("Optimization Cycle:", iteration+1)						# For developers, iteration starts at 0, for users iteration starts at 1 (so print that)
		normal_H = hamiltonian(Enuc, h1, V2, fluctuations)
		print(normal_H._C[0])
		Hop  = Cspace.lin_op(super_operator_wrapper(normal_H, resources))
		Heff = Pop | Hop
		textlog(indent("Normal-ordered H built and projected into {}CIRS space.".format("rotated " if iteration>0 else "")))
		Cvec = Cspace.member(operator(fluctuations,C=1))					# since only constant part is 1, this effectively represents the reference
		textlog(indent("Lanczos guess set to {}tensor-product reference state.".format(rotate_scf.antisymmetrized)))
		#
		textlog(indent("Extracting lowest-energy CIRS state ..."))
		evalue,evec = lanczos.lowest_eigen(Heff, Cvec, thresh=thresh, dim=krylov_dim)
		excitations = excite_mask(evec.v)
		rot_norm = excitations.dot(excitations)
		textlog(indent("Norm of singles projection of lowest state =", rot_norm))
		#
		textlog(indent("Using singles amplitudes to rotate one-body reference functions ..."))
		E0, E1, energy_diffs, fluctuations, h1, V2 = rotate_scf(h1, V2, excitations)
		textlog(indent("-------------------------"))
		textlog(indent("E_SCF = {}   (sum of one-body energies:  E_0 = {})".format(E1,E0)))
		textlog(indent("-------------------------"))
		#
		iteration += 1
	textlog("==================================================================================")
	textlog("Generalized SCF converged! ({} cycles)".format(iteration))
	textlog("==================================================================================")
	normal_H = hamiltonian(Enuc, h1, V2, fluctuations)
	textlog("Unprojected normal-ordered H built in final CIRS space.")
	textlog("==================================================================================")
	return fluctuations, normal_H, energy_diffs

def ccsd(fluctuations, H, invF, textlog, resources):
	out = output(log=textlog)
	#
	T = nested_operator(fluctuations)		# Initialized to zero
	out.log("T initialized.")
	BCH = BakerCampbellHausdorff(fluctuations, resources)
	out.log("Baker-Campbell-Hausdorff module initialized.")
	#
	out.log("Performing coupled-cluster iterations ...")
	t_ccsd_start = time.time()
	E = CCSD(H, T, invF, BCH, space_traits, out.log.sublog())
	t_ccsd_end = time.time()
	#
	out.log("CCSD Energy =", E, 'Hartree')
	out.log("CCSD TIME {}".format(t_ccsd_end - t_ccsd_start))
	return out(energy=E)



### Helper function

def parse_input():
	E_nuc, num_alpha_elec, num_beta_elec, num_spin_orb = 0., 1, 1, 8
	#
	T_mat, N_mat, V_mat_raw = scf_routine_hacked.main("1.0.in")
	h_mat = T_mat + N_mat
	V_mat = digest_V_mat(V_mat_raw, num_spin_orb)
	#
	occ_orbs = []
	vrt_orbs = []
	for p in range(num_alpha_elec):  occ_orbs += [p]
	for p in range(num_beta_elec):   occ_orbs += [p+(num_spin_orb//2)]
	for p in range(num_spin_orb):
		if p not in occ_orbs:  vrt_orbs += [p]
	#
	return output(Enuc=E_nuc, h1=h_mat, V2=V_mat), output(occ=occ_orbs, vrt=vrt_orbs)



#######
main()# Do it!
#######
