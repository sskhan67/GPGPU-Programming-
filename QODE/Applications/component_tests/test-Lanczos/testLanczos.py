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
import time
import numpy
import qode



def io_setup(argv):
    debug = False
    numpy.set_printoptions(linewidth=200)
    data = qode.util.read_input.from_file(argv[0]+'/params.py')
    data = qode.util.read_input.from_argv(argv[1:],write_to=data)
    mat  = qode.util.random_matrix(data.dim)
    vecs = [ qode.util.basis_vector(i,data.dim) for i in data.guess ]
    if debug:
        print(mat)
        print(vecs)
    return mat,vecs,data.thresh,data.psize,data.block_action

def control(mat,n):
    debug = False
    eigenvals,eigenvecs = numpy.linalg.eigh(mat)
    if debug:
        print(eigenvals)
        diag = eigenvecs.T * mat * eigenvecs
        print(diag)
    print("Lowest eigenvalues from full diagonalization:  {}".format(sorted(eigenvals)[:n]))

def experiment(mat,vecs,thresh,psize,block_action):
    traits = qode.math.numpy_space.real_traits(len(vecs[0]))
    S = qode.math.linear_inner_product_space(traits)
    H = S.lin_op(mat)
    v = [ S.member(vec) for vec in vecs ]
    eigenpairs = qode.math.lanczos.lowest_eigen(H,v,thresh,psize,block_action=block_action)
    print("Lowest eigenvalues from Lanczos extraction:    {}".format(list(zip(*eigenpairs))[0]))



def main(argv):
    debug = False
    t0 = time.time()
    #
    mat,vecs,thresh,psize,block_action = io_setup(argv)
    t1 = time.time()
    #
    if debug:
        print("Skipping full LAPACK diagonalization.")
        t2 = t1
    else:
        control(mat,len(vecs))
        t2 = time.time()
    #
    experiment(mat,vecs,thresh,psize,block_action)
    t3 = time.time()
    #
    tmat,tdiag,tlanc,ttot = (t1-t0)/60 , (t2-t1)/60 , (t3-t2)/60 , (t3-t0)/60
    print("matrix build:          {} min".format(tmat))
    print("full diagonalization:  {} min".format(tdiag))
    print("lanczos extraction:    {} min".format(tlanc))
    print("Total:                 {} min".format(ttot))

if __name__ ==  "__main__":  main(sys.argv[1:])
