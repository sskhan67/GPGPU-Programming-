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
import numpy
from   numerical_engines     import H_contractions as contract	# The functions in here are the present optimization targets
from   build_density_tensors import build_density_tensors



# Right now, most of the code outside this file is set up for only having two fragments.  Since many of the
# calculations here are so repetitive in structure for larger systems, we can mimic a full-size calculation at some
# places, by just doing the same calculation many times (SIMD on the GPUs).  The size parameter n_frag here
# specifies the number of such "dummy" fragments and really only changes some loop lengths below.
#############
n_frag = 10
#############



# Presently, we are testing only on homogeneous systems (all fragments are chemically identical).
# Here are some (temporarily) hard-coded system parameters for frozen-core 6-31G Be atoms.
n_elec_frag = 4    # The number of electrons when the fragment is neutral
n_orb_frag  = 9    # The number of *spatial* orbitals for each fragment
n_core_frag = 1    # The number of frozen (uncorrelated) core orbitals



# Here we load up the system-wide primitive description of the physical interactions of individual electrons.
# These will be contracted with descriptions of the fragment states below to give a compressed ("renormalized")
# description of interactions between entire fragments.  This is the point at which we are really limited
# by external code constraints.  We are presently only generating these files for dimers since we are not being
# very clever about how the tensors are stored.  It is also presently hardcoded to look for specific systems
# (Be dimers) only varying the distance between them.
h = numpy.load("data/semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_h.npy".format(sys.argv[1]))    # 1-electron integrals (kinetic energy, nuclear attraction)
V = numpy.load("data/semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_V.npy".format(sys.argv[1]))    # 2-electron integrals (electron repulsion)
V = (1/4.) * ( V - numpy.swapaxes(V,2,3) )                                                        # Account for electron indistinguishability at integrals level
print("integrals loaded")



# Load descriptions of the fragment many-electron states.  Then process them (using build_density_tensors) to arrive at
# quantities suitable for contracting with the integrals loaded above.  Since we are presently only working with
# homogeneous systems, only one "example" fragment is loaded.  We can/will extend the code up to quite large (and varied)
# homogeneous systems, without having to change this step.
class states(object): pass      #
neutral = states()              # Small dictionaries to hold symbolic representations of fragment electron configurations
cation  = states()              # and their associated numerical coefficients in the fragment wavefunctions.
anion   = states()              #
neutral.coeffs  = numpy.load("data/atomic_states/H/16-115-550/{}/Z_2e.npy".format(sys.argv[2])).T   #
neutral.configs = numpy.load("data/atomic_states/data/configs_2e-stable.npy")                       # In a supersystem, any given
cation.coeffs   = numpy.load("data/atomic_states/H/16-115-550/{}/Z_1e.npy".format(sys.argv[2])).T   # fragment may be found as neutral
cation.configs  = numpy.load("data/atomic_states/data/configs_1e-stable.npy")                       # or with an extra or missing electron
anion.coeffs    = numpy.load("data/atomic_states/H/16-115-550/{}/Z_3e.npy".format(sys.argv[2])).T   # (anionic or cationic, respectively).
anion.configs   = numpy.load("data/atomic_states/data/configs_3e-stable.npy")                       #
rho_frag = build_density_tensors({0:neutral, +1:cation, -1:anion}, n_orb_frag, n_core_frag)   # Compute the density tensors
charge_states_frag = { 0:neutral.coeffs.shape[0], +1:cation.coeffs.shape[0], -1:anion.coeffs.shape[0]}   # Number of fragment states for each fragment charge (0,+1,-1)
print("rho computed")    # "rho" is the symbol traditionally given to a density matrix or tensor



# This class manages the many different tensor contractions of primitive integrals and fragment density tensors (loaded above),
# which need to be performed in order to build the renormalized supersystem Hamiltonian matrix.  Its member functions are called
# from inside loops over fragments and states of those fragments (below).  This presently builds matrix elements one at a
# time.  Since the building of each matrix element is structurally similar (just a call to one of the member functions
# here, with information utlimatley pointing to different tensors, or blocks thereof), there should be good SIMD potential.
class build_matrix_elements(object):
    def __init__(self, rho, h, V, n_orb, n_elec):    # n_elec refers to the fragments in their reference states
        n_orb = [2*n for n in n_orb]                 # convert to spin orbitals (does not alter original list)
        n_orb_tot = sum(n_orb)
        n_orb_tot = sum(n_orb[:2])                   # Nasty hack until this really does more than just dimers
        self.data = rho, h, V, n_orb, n_orb_tot, n_elec
    def monomer(self, fragment, I, J):
        rho, h, V, n_orb, n_orb_tot, n_elec = self.data
        m = fragment
        i_m_chg, i_m_idx = I
        j_m_chg, j_m_idx = J
        if i_m_chg == j_m_chg:
            Rca    = rho[m]["ca"  ][(i_m_chg, j_m_chg)][i_m_idx][j_m_idx]
            Rccaa  = rho[m]["ccaa"][(i_m_chg, j_m_chg)][i_m_idx][j_m_idx]
            offset = sum(n_orb[:m])
            n_orb  =     n_orb[m]
            return contract.monomer(offset, n_orb_tot, n_orb, Rca, Rccaa, h, V)
        else:
            return 0
    def dimer(self, fragments, I, J):
        rho, h, V, n_orb, n_orb_tot, n_elec = self.data
        m1, m2 = fragments
        i_m1, i_m2 = I
        j_m1, j_m2 = J
        i_m1_chg, i_m1_idx = i_m1
        i_m2_chg, i_m2_idx = i_m2
        j_m1_chg, j_m1_idx = j_m1
        j_m2_chg, j_m2_idx = j_m2
        if (i_m1_chg+i_m2_chg == j_m1_chg+j_m2_chg) and (i_m1_chg >= j_m1_chg-2) and (i_m1_chg <= j_m1_chg+2):
            offset1 = sum(n_orb[:m1])
            offset2 = sum(n_orb[:m2])
            n_orb1  =     n_orb[m1]
            n_orb2  =     n_orb[m2]
            if i_m1_chg == j_m1_chg-2:
                Rcc1 = rho[m1]["cc"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Raa2 = rho[m2]["aa"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                return contract.dimer_2min2pls(offset1, offset2, n_orb_tot, n_orb1, n_orb2, Rcc1, Raa2, h, V)
            if i_m1_chg == j_m1_chg-1:
                Rc1   = rho[m1]["c"  ][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rcca1 = rho[m1]["cca"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Ra2   = rho[m2]["a"  ][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                Rcaa2 = rho[m2]["caa"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                sign = -(-1)**((n_elec[m2]+j_m2_chg) % 2)
                return sign * contract.dimer_1min1pls(offset1, offset2, n_orb_tot, n_orb1, n_orb2, Rc1, Rcca1, Ra2, Rcaa2, h, V)
            if i_m1_chg == j_m1_chg:
                Rca1 = rho[m1]["ca"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rca2 = rho[m2]["ca"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                return contract.dimer_00(offset1, offset2, n_orb_tot, n_orb1, n_orb2, Rca1, Rca2, h, V)
            if i_m1_chg == j_m1_chg+1:
                Ra1   = rho[m1]["a"  ][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rcaa1 = rho[m1]["caa"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rc2   = rho[m2]["c"  ][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                Rcca2 = rho[m2]["cca"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                sign = +(-1)**((n_elec[m2]+j_m2_chg) % 2)
                return sign * contract.dimer_1pls1min(offset1, offset2, n_orb_tot, n_orb1, n_orb2, Ra1, Rcaa1, Rc2, Rcca2, h, V)
            if i_m1_chg == j_m1_chg+2:
                Raa1 = rho[m1]["aa"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]
                Rcc2 = rho[m2]["cc"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]
                return contract.dimer_2pls2min(offset1, offset2, n_orb_tot, n_orb1, n_orb2, Raa1, Rcc2, h, V)
        else:
                return 0



# Here we set up the homogeneous supersystem by simply making many pointers to the information for a single fragment.
# The information that distinguishes among them (e.g., who is next to to whom) is carried by the interaction
# integrals.  While we will need the actual integrals to do large systems, this setup will not change until we
# start doing inhomogeneous systems.
n_elec        =        [n_elec_frag]*n_frag
n_orb         =         [n_orb_frag]*n_frag
rho           =           [rho_frag]*n_frag
charge_states = [charge_states_frag]*n_frag



# Instantiate the engine that computes the matrix elements.  n_orb will eventually not be necessary,
# when the integrals get passed in a more structured way that exposes internal blocks.
compute = build_matrix_elements(rho, h, V, n_orb, n_elec)



# Allocate storage for all of the matrix elements to be computed
frag_dims = [sum(n for chg,n in chgstM.items()) for chgstM in charge_states]
H1 = [numpy.zeros((dim,dim)) for dim in frag_dims]
H2 = [[numpy.zeros((dimM*dimM, dimN*dimN)) for dimN in frag_dims[:M]] for M,dimM in enumerate(frag_dims)]

# Compute the renormalized monomer terms
for M in range(n_frag):
    i = 0
    for bra_chg,bra_states in charge_states[M].items():
        for bra in range(bra_states):
            I = bra_chg, bra
            j = 0
            for ket_chg,ket_states in charge_states[M].items():
                for ket in range(ket_states):
                    J = ket_chg, ket
                    #H1[M][i,j] = compute.monomer(M, I, J)    # This would presently crash because there are only enough integrals for two systems (0 and 1) ...
                    H1[M][i,j] = compute.monomer(0, I, J)     # ... but this will give the same answer, since the system is presently homogeneous.
                    j += 1
            i += 1
    print("Finished H1[{}]".format(M))

# Compute the renormalized dimer terms
for M in range(n_frag):
    for N in range(M):
        i = 0
        for bra_chgM,bra_statesM in charge_states[M].items():
            for braM in range(bra_statesM):
                for bra_chgN,bra_statesN in charge_states[N].items():
                    for braN in range(bra_statesN):
                        I = (bra_chgM,braM), (bra_chgN,braN)
                        j = 0
                        for ket_chgM,ket_statesM in charge_states[M].items():
                            for ketM in range(ket_statesM):
                                for ket_chgN,ket_statesN in charge_states[N].items():
                                    for ketN in range(ket_statesN):
                                        J = (ket_chgM,ketM), (ket_chgN,ketN)
                                        #H2[M][N][i,j] = compute.dimer((M,N), I, J)    # This would presently crash because there are only enough integrals for the 0-1 dimer ...
                                        H2[M][N][i,j] = compute.dimer((0,1), I, J)     # ... but this will do a calculation of the same structure, just with the wrong integrals.
                                        j += 1
                        i += 1
        print("Finished H2[{}][{}]".format(M,N))

# Dump the renormalized Hamiltonian to disk for use in XR-CC calculation
#numpy.save("H1.npy", H1)
#numpy.save("H2.npy", H2)
