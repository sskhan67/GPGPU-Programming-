#    (C) Copyright 2019 Anthony D. Dutoi
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

# Use within a Psi4 conda environment

import numpy
import psi4
from   qode.util.PyC import Double
import qode.atoms.integrals.spatial_to_spin            as spatial_to_spin
import qode.atoms.integrals.external_engines.psi4_ints as integrals

# Hard-coded assumptions
basis = "6-31G"
system = "Be"	# homodimer



def biorthogonalize(H, V, S):
    S_inv = numpy.linalg.inv(S) 
    I = S_inv @ S
    refine = numpy.linalg.inv(I)
    S_inv = refine @ S_inv
    H = S_inv @ H
    for _ in range(2):  V = numpy.tensordot(S_inv, V, axes=([1],[1]))	# cycle through first two tensor axes
    return H, V

def to_semiMO_basis(S, H, V, C1):
    n1 = C1.shape[0]
    N  =  S.shape[0] // n1
    C  = numpy.zeros((N*n1, N*n1))
    for i in range(0, N*n1, n1):  C[i:i+n1, i:i+n1] = C1
    S = C.T @ S @ C
    H = C.T @ H @ C
    for _ in range(4):  V = numpy.tensordot(V, C, axes=([0],[0]))	# cycle through the tensor axes (this assumes everything is real)
    return S, H, V




def block_ints(dist_str, Cmat1):
    sys1 = "{sys}\n".format(sys=system)
    S1, T1, U1, V1, _ = integrals.AO_ints(sys1, basis)
    n_orb_frag = S1.shape[0]

    sys2 = "{sys}\n{sys} 1 {dist}\n".format(sys=system,dist=dist_str)
    S2, T2, U2, V2, _ = integrals.AO_ints(sys2, basis)
    H2 = T2 + U2

    # Biorthogonal semi-MO SCF Be chain
    S_semiMO, H_semiMO, V_semiMO = to_semiMO_basis(S2, H2, V2, Cmat1)
    H_semiMO_bi, V_semiMO_bi = biorthogonalize(H_semiMO, V_semiMO, S_semiMO)

    # Integrals blocked by fragment, then spin (alpha block before beta block)
    H_tmp, V_tmp = spatial_to_spin.one_electron(H_semiMO_bi), spatial_to_spin.two_electron(V_semiMO_bi)
    if H_tmp.shape[0] != 4*n_orb_frag:  raise LogicError
    H_semiMO_bi = numpy.zeros((4*n_orb_frag,4*n_orb_frag))
    for p1,p2 in [(0,0),(n_orb_frag,2*n_orb_frag),(2*n_orb_frag,n_orb_frag),(3*n_orb_frag,3*n_orb_frag)]:
        for p in range(n_orb_frag):
            for q1,q2 in [(0,0),(n_orb_frag,2*n_orb_frag),(2*n_orb_frag,n_orb_frag),(3*n_orb_frag,3*n_orb_frag)]:
                for q in range(n_orb_frag):
                    H_semiMO_bi[p1+p, q1+q] = H_tmp[p2+p, q2+q]
    V_semiMO_bi = numpy.zeros((4*n_orb_frag,4*n_orb_frag,4*n_orb_frag,4*n_orb_frag))
    for p1,p2 in [(0,0),(n_orb_frag,2*n_orb_frag),(2*n_orb_frag,n_orb_frag),(3*n_orb_frag,3*n_orb_frag)]:
        for p in range(n_orb_frag):
            for q1,q2 in [(0,0),(n_orb_frag,2*n_orb_frag),(2*n_orb_frag,n_orb_frag),(3*n_orb_frag,3*n_orb_frag)]:
                for q in range(n_orb_frag):
                    for r1,r2 in [(0,0),(n_orb_frag,2*n_orb_frag),(2*n_orb_frag,n_orb_frag),(3*n_orb_frag,3*n_orb_frag)]:
                        for r in range(n_orb_frag):
                            for s1,s2 in [(0,0),(n_orb_frag,2*n_orb_frag),(2*n_orb_frag,n_orb_frag),(3*n_orb_frag,3*n_orb_frag)]:
                                for s in range(n_orb_frag):
                                    V_semiMO_bi[p1+p, q1+q, r1+r, s1+s] = V_tmp[p2+p, q2+q, r2+r, s2+s]
                                                  
    # Account for electron indistinguishability at integrals level
    V_semiMO_bi = (1/4.) * ( V_semiMO_bi - numpy.swapaxes(V_semiMO_bi,2,3) )

    N = 2*n_orb_frag
    h = [[None,None],[None,None]]
    h[0][0] = numpy.array(H_semiMO_bi[:N,:N], dtype=Double.numpy)
    h[0][1] = numpy.array(H_semiMO_bi[:N,N:], dtype=Double.numpy)
    h[1][0] = numpy.array(H_semiMO_bi[N:,:N], dtype=Double.numpy)
    h[1][1] = numpy.array(H_semiMO_bi[N:,N:], dtype=Double.numpy)
    V = [[[[None,None],[None,None]],[[None,None],[None,None]]],[[[None,None],[None,None]],[[None,None],[None,None]]]]
    V[0][0][0][0] = numpy.array(V_semiMO_bi[:N,:N,:N,:N], dtype=Double.numpy)
    V[0][0][0][1] = numpy.array(V_semiMO_bi[:N,:N,:N,N:], dtype=Double.numpy)
    V[0][0][1][0] = numpy.array(V_semiMO_bi[:N,:N,N:,:N], dtype=Double.numpy)
    V[0][0][1][1] = numpy.array(V_semiMO_bi[:N,:N,N:,N:], dtype=Double.numpy)
    V[0][1][0][0] = numpy.array(V_semiMO_bi[:N,N:,:N,:N], dtype=Double.numpy)
    V[0][1][0][1] = numpy.array(V_semiMO_bi[:N,N:,:N,N:], dtype=Double.numpy)
    V[0][1][1][0] = numpy.array(V_semiMO_bi[:N,N:,N:,:N], dtype=Double.numpy)
    V[0][1][1][1] = numpy.array(V_semiMO_bi[:N,N:,N:,N:], dtype=Double.numpy)
    V[1][0][0][0] = numpy.array(V_semiMO_bi[N:,:N,:N,:N], dtype=Double.numpy)
    V[1][0][0][1] = numpy.array(V_semiMO_bi[N:,:N,:N,N:], dtype=Double.numpy)
    V[1][0][1][0] = numpy.array(V_semiMO_bi[N:,:N,N:,:N], dtype=Double.numpy)
    V[1][0][1][1] = numpy.array(V_semiMO_bi[N:,:N,N:,N:], dtype=Double.numpy)
    V[1][1][0][0] = numpy.array(V_semiMO_bi[N:,N:,:N,:N], dtype=Double.numpy)
    V[1][1][0][1] = numpy.array(V_semiMO_bi[N:,N:,:N,N:], dtype=Double.numpy)
    V[1][1][1][0] = numpy.array(V_semiMO_bi[N:,N:,N:,:N], dtype=Double.numpy)
    V[1][1][1][1] = numpy.array(V_semiMO_bi[N:,N:,N:,N:], dtype=Double.numpy)

    return h, V
