#    (C) Copyright 2018, 2019 Anthony D. Dutoi
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
import pickle
import numpy
from numerical_engines import H_contractions as contract
from build_density_tensors import build_density_tensors



h = numpy.load("data/semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_h.npy".format(sys.argv[1]))	# Pray that the system . . .
V = numpy.load("data/semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_V.npy".format(sys.argv[1]))	# . . . is in the correct orientation.
V = (1/4.) * ( V - numpy.swapaxes(V,2,3) )
E_nuc = 4 * 4 * 0.52917721067 / float(sys.argv[1])

print("integrals loaded")


class states(object): pass

neutral = states()
neutral.coeffs  = numpy.load("data/atomic_states/H/16-115-550/{}/Z_2e.npy".format(sys.argv[2])).T
neutral.configs = numpy.load("data/atomic_states/data/configs_2e-stable.npy")
cation = states()
cation.coeffs  = numpy.load("data/atomic_states/H/16-115-550/{}/Z_1e.npy".format(sys.argv[2])).T
cation.configs = numpy.load("data/atomic_states/data/configs_1e-stable.npy")
anion = states()
anion.coeffs  = numpy.load("data/atomic_states/H/16-115-550/{}/Z_3e.npy".format(sys.argv[2])).T
anion.configs = numpy.load("data/atomic_states/data/configs_3e-stable.npy")

z_lists = { 0:neutral , +1:cation, -1:anion }
rho_atom = build_density_tensors(z_lists, n_orbs=9, n_core=1)	# n_orbs and n_core are numbers of spatial orbitals
rho = [rho_atom, rho_atom]	# homo-dimer

print("rho computed")

class build_matrix_elements(object):
    def __init__(self, rho, h, V, n_orb_sys, n_elec_sys):	# n_elec_sys refers to the fragments in their reference states
        n_orb_tot = sum(n_orb_sys)
        self.data = rho, h, V, n_orb_sys, n_orb_tot, n_elec_sys
    def monomer(self, subsystem, I, J):
        rho, h, V, n_orb_sys, n_orb_tot, n_elec_sys = self.data
        m = subsystem
        i_m_chg, i_m_idx = I
        j_m_chg, j_m_idx = J
        if i_m_chg == j_m_chg:
            Rca       = rho[m]["ca"  ][(i_m_chg, j_m_chg)][i_m_idx][j_m_idx]
            Rccaa     = rho[m]["ccaa"][(i_m_chg, j_m_chg)][i_m_idx][j_m_idx]
            offset    = sum(n_orb_sys[:m])
            n_orb_sys =     n_orb_sys[m]
            return contract.monomer(offset, n_orb_tot, n_orb_sys, Rca, Rccaa, h, V)
        else:
            return 0
    def dimer(self, subsystems, I, J):
        rho, h, V, n_orb_sys, n_orb_tot, n_elec_sys = self.data
        m1, m2 = subsystems
        i_m1, i_m2 = I
        j_m1, j_m2 = J
        i_m1_chg, i_m1_idx = i_m1
        i_m2_chg, i_m2_idx = i_m2
        j_m1_chg, j_m1_idx = j_m1
        j_m2_chg, j_m2_idx = j_m2
        if (i_m1_chg+i_m2_chg == j_m1_chg+j_m2_chg) and (i_m1_chg >= j_m1_chg-2) and (i_m1_chg <= j_m1_chg+2):
            offset1    = sum(n_orb_sys[:m1])
            offset2    = sum(n_orb_sys[:m2])
            n_orb_sys1 =     n_orb_sys[m1]
            n_orb_sys2 =     n_orb_sys[m2]
            if i_m1_chg == j_m1_chg-2:
                Rcc1 = rho[m1]["cc"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Raa2 = rho[m2]["aa"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                return contract.dimer_2min2pls(offset1, offset2, n_orb_tot, n_orb_sys1, n_orb_sys2, Rcc1, Raa2, h, V)
            if i_m1_chg == j_m1_chg-1:
                Rc1   = rho[m1]["c"  ][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rcca1 = rho[m1]["cca"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Ra2   = rho[m2]["a"  ][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                Rcaa2 = rho[m2]["caa"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                sign = -(-1)**((n_elec_sys[m2]+j_m2_chg) % 2)
                return sign * contract.dimer_1min1pls(offset1, offset2, n_orb_tot, n_orb_sys1, n_orb_sys2, Rc1, Rcca1, Ra2, Rcaa2, h, V)
            if i_m1_chg == j_m1_chg:
                Rca1 = rho[m1]["ca"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rca2 = rho[m2]["ca"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                return contract.dimer_00(offset1, offset2, n_orb_tot, n_orb_sys1, n_orb_sys2, Rca1, Rca2, h, V)
            if i_m1_chg == j_m1_chg+1:
                Ra1   = rho[m1]["a"  ][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rcaa1 = rho[m1]["caa"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rc2   = rho[m2]["c"  ][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                Rcca2 = rho[m2]["cca"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                sign = +(-1)**((n_elec_sys[m2]+j_m2_chg) % 2)
                return sign * contract.dimer_1pls1min(offset1, offset2, n_orb_tot, n_orb_sys1, n_orb_sys2, Ra1, Rcaa1, Rc2, Rcca2, h, V)
            if i_m1_chg == j_m1_chg+2:
                Raa1 = rho[m1]["aa"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rcc2 = rho[m2]["cc"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                return contract.dimer_2pls2min(offset1, offset2, n_orb_tot, n_orb_sys1, n_orb_sys2, Raa1, Rcc2, h, V)
        else:
                return 0


n_orb_tot = h.shape[0]

compute = build_matrix_elements(rho, h, V, [n_orb_tot//2, n_orb_tot//2], [4,4])

frag_dim = sum([z_lists[chg].coeffs.shape[0] for chg in (0,+1,-1)])

H1_0 = numpy.zeros((frag_dim, frag_dim))
i = 0
for bra_chg0 in (0,+1,-1):
    bra_n0 = z_lists[bra_chg0].coeffs.shape[0]
    for bra0 in range(bra_n0):
        I = bra_chg0, bra0
        j = 0
        for ket_chg0 in (0,+1,-1):
            ket_n0 = z_lists[ket_chg0].coeffs.shape[0]
            for ket0 in range(ket_n0):
                J = ket_chg0, ket0
                H1_0[i,j] = compute.monomer(0, I, J)
                j += 1
        i += 1

H1_1 = numpy.zeros((frag_dim, frag_dim))
i = 0
for bra_chg1 in (0,+1,-1):
    bra_n1 = z_lists[bra_chg1].coeffs.shape[0]
    for bra1 in range(bra_n1):
        I = bra_chg1, bra1
        j = 0
        for ket_chg1 in (0,+1,-1):
            ket_n1 = z_lists[ket_chg1].coeffs.shape[0]
            for ket1 in range(ket_n1):
                J = ket_chg1, ket1
                H1_1[i,j] = compute.monomer(1, I, J)
                j += 1
        i += 1

H2_01   = numpy.zeros((frag_dim*frag_dim, frag_dim*frag_dim))
H_dimer = numpy.zeros((frag_dim*frag_dim, frag_dim*frag_dim))
i = 0
for bra_chg0 in (0,+1,-1):
    bra_n0 = z_lists[bra_chg0].coeffs.shape[0]
    for bra0 in range(bra_n0):
        for bra_chg1 in (0,+1,-1):
            bra_n1 = z_lists[bra_chg1].coeffs.shape[0]
            for bra1 in range(bra_n1):
                I = (bra_chg0,bra0), (bra_chg1,bra1)
                j = 0
                for ket_chg0 in (0,+1,-1):
                    ket_n0 = z_lists[ket_chg0].coeffs.shape[0]
                    for ket0 in range(ket_n0):
                        for ket_chg1 in (0,+1,-1):
                            ket_n1 = z_lists[ket_chg1].coeffs.shape[0]
                            for ket1 in range(ket_n1):
                                J = (ket_chg0,ket0), (ket_chg1,ket1)
                                coupling = compute.dimer((0,1), I, J)
                                H2_01[i,j]   = coupling
                                H_dimer[i,j] = coupling
                                if I[0]==J[0]:  H_dimer[i,j] += compute.monomer(1, I[1], J[1])
                                if I[1]==J[1]:  H_dimer[i,j] += compute.monomer(0, I[0], J[0])
                                j += 1
                i += 1

# did this at the very last minute
#numpy.save("H1_0.npy",  H1_0)
#numpy.save("H1_1.npy",  H1_1)
#numpy.save("H2_01.npy", H2_01)






# Extract the lowest eigenvalue from the compressed Hamiltonian for comparison
print("Diagonalizing ... ", end="", flush=True)
Evals,Evecs = numpy.linalg.eig(H_dimer)
Evals,Evecs = zip(*sorted(zip( Evals, Evecs.T.tolist() ), key=lambda p: p[0]))
print("Done.  \n\nE_gs  = {}  [compressed H]\n".format(Evals[0] + E_nuc))

"""\
idx = 0
for ket_chg0 in (0,+1,-1):
    ket_n0 = z_lists[ket_chg0].coeffs.shape[0]
    for ket0 in range(ket_n0):
        for ket_chg1 in (0,+1,-1):
            ket_n1 = z_lists[ket_chg1].coeffs.shape[0]
            for ket1 in range(ket_n1):
                val = " "
                if Evecs[0][idx]>1e-2:  val = "{: f}".format(Evecs[0][idx])
                print("{:2d} {:2d} {:2d} {:2d} {:s}".format(ket_chg0, ket0, ket_chg1, ket1, val))
                idx += 1
"""



"""\
# Check our work
out, resources  = output(log=textlog(echo=True)), parallel.resources(1)
out.log("EXCITONIC CCSD CALCULATION")
filestem_pattern = path+"/{0}/{1}"
excitonic.ccsd(hamiltonian.load_subH(filestem_pattern, [[".",dist],[".","."]]), out, resources)

f_mat = numpy.load(path+"/"+dist+"/f.npy")     # Could call fci.main() here, but until everything is rock solid, I like independence of this simple dimer-only H-build
V_mat = numpy.load(path+"/"+dist+"/V.npy")
out.log("Building FCI matrix ...")
for i in range(states_per_frag):
	for j in range(states_per_frag):
		Mij = i*states_per_frag + j
		for k in range(states_per_frag):
			for l in range(states_per_frag):
				Nkl = k*states_per_frag + l
				if i==k:  V_mat[Mij][Nkl] += f_mat[j][l]
				if j==l:  V_mat[Mij][Nkl] += f_mat[i][k]
out.log("Diagonalizing ...")
vals, vecs = numpy.linalg.eig(V_mat)
vals, vecs = zip(*sorted(zip( vals, vecs.T.tolist() ), key=lambda p: p[0]))
out.log("FCI Energy  =", vals[0])

out.log.write(path+"/"+dist+"/log.txt")
"""
