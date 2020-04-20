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
import numpy
import qode
from   qode.math import sqrt
from   qode.math import linear_inner_product_space as space
from   qode.math import numpy_space                as np_space



def main(argv):

    data   = qode.util.read_input.from_file(argv[0]+'/params.py',namespace=globals())
    basisA = qode.util.basis_vector(data.basis[0],data.dim)
    basisB = qode.util.basis_vector(data.basis[1],data.dim)
    basisC = qode.util.basis_vector(data.basis[2],data.dim)
    vec1 = data.vec1[0]*basisA + data.vec1[1]*basisB
    vec2 = data.vec2[0]*basisA + data.vec2[1]*basisB
    print(vec1)
    print(vec2)

    S = space( np_space.real_traits(data.dim) )
    vecs = qode.math.vector_set(S)

    basisA = S.member(basisA)
    basisB = S.member(basisB)
    basisC = S.member(basisC)
    vec1   = S.member(vec1)
    vec2   = S.member(vec2)

    vecs.append(vec1)
    vecs.append(vec2)
    print(vecs.overlaps())
    print(vecs.is_orthonormal())

    print(vecs.projections(basisA))
    print(vecs.projections(basisB))

    print(vecs.project(basisA))
    print(vecs.project(basisB))
    print(vecs.project(basisC))

    P = vecs.projection_opr()
    print( (P|basisA).v )
    print( (P|basisB).v )
    print( (P|basisC).v )

    vecs.append(basisA)
    vecs.append(basisB)
    print(vecs.overlaps())
    print(vecs.is_orthonormal())



if __name__ ==  "__main__":  main(sys.argv[1:])


