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

import numpy
from qode.util.PyC import Double

def block_ints(dist_string, n_orb_frag):

    hraw = numpy.load("semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_h.npy".format(dist_string))    # 1-electron integrals (kinetic energy, nuclear attraction)
    Vraw = numpy.load("semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_V.npy".format(dist_string))    # 2-electron integrals (electron repulsion)
    Vraw = (1/4.) * ( Vraw - numpy.swapaxes(Vraw,2,3) )                                                  # Account for electron indistinguishability at integrals level

    N = 2*n_orb_frag
    h = [[None,None],[None,None]]
    h[0][0] = numpy.array(hraw[:N,:N], dtype=Double.numpy)
    h[0][1] = numpy.array(hraw[:N,N:], dtype=Double.numpy)
    h[1][0] = numpy.array(hraw[N:,:N], dtype=Double.numpy)
    h[1][1] = numpy.array(hraw[N:,N:], dtype=Double.numpy)
    V = [[[[None,None],[None,None]],[[None,None],[None,None]]],[[[None,None],[None,None]],[[None,None],[None,None]]]]
    V[0][0][0][0] = numpy.array(Vraw[:N,:N,:N,:N], dtype=Double.numpy)
    V[0][0][0][1] = numpy.array(Vraw[:N,:N,:N,N:], dtype=Double.numpy)
    V[0][0][1][0] = numpy.array(Vraw[:N,:N,N:,:N], dtype=Double.numpy)
    V[0][0][1][1] = numpy.array(Vraw[:N,:N,N:,N:], dtype=Double.numpy)
    V[0][1][0][0] = numpy.array(Vraw[:N,N:,:N,:N], dtype=Double.numpy)
    V[0][1][0][1] = numpy.array(Vraw[:N,N:,:N,N:], dtype=Double.numpy)
    V[0][1][1][0] = numpy.array(Vraw[:N,N:,N:,:N], dtype=Double.numpy)
    V[0][1][1][1] = numpy.array(Vraw[:N,N:,N:,N:], dtype=Double.numpy)
    V[1][0][0][0] = numpy.array(Vraw[N:,:N,:N,:N], dtype=Double.numpy)
    V[1][0][0][1] = numpy.array(Vraw[N:,:N,:N,N:], dtype=Double.numpy)
    V[1][0][1][0] = numpy.array(Vraw[N:,:N,N:,:N], dtype=Double.numpy)
    V[1][0][1][1] = numpy.array(Vraw[N:,:N,N:,N:], dtype=Double.numpy)
    V[1][1][0][0] = numpy.array(Vraw[N:,N:,:N,:N], dtype=Double.numpy)
    V[1][1][0][1] = numpy.array(Vraw[N:,N:,:N,N:], dtype=Double.numpy)
    V[1][1][1][0] = numpy.array(Vraw[N:,N:,N:,:N], dtype=Double.numpy)
    V[1][1][1][1] = numpy.array(Vraw[N:,N:,N:,N:], dtype=Double.numpy)

    return h, V
