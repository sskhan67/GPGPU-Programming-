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
    vec  = qode.util.column_matrix(data.guess,data.dim)
    if debug:
        print(mat)
        print(vec)
    return mat,vec,data.thresh,data.psize

def control(mat):
    debug = False
    eigenvals,eigenvecs = numpy.linalg.eigh(mat)
    lowest = qode.util.min.index(eigenvals)
    evalue = eigenvals[lowest]
    print("Lowest eigenvalue from full diagonalization:  {}".format(evalue))
    if debug:
        print(lowest)
        diag = eigenvecs.T * mat * eigenvecs
        print(diag)

def experiment(mat,vec,thresh,psize):
    H     = qode.math.numpy_space.operator(mat)
    print(H)
    print(H.mat)
    guess = qode.math.numpy_space.vector(vec)
    evalue,evec = qode.math.lanczos.lowest_eigen(H,guess,thresh,psize)
    print("Lowest eigenvalue from Lanczos extraction:    {}".format(evalue))



def main(argv):
    debug = True
    t0 = time.time()
    #
    # mat,vec,thresh,psize = io_setup(argv)
    mat = numpy.matrix( [[1., 0.1, 0.1, 0.1], [0.1,2.0,0.1,0.1],[0.1,0.1,3.0,0.1],[0.1,0.1,0.1,4.0]] )
    vec = numpy.matrix( [[1.0],[0.0],[0.0],[0.0]] )
    thresh = 1e-6
    psize = 4
    t1 = time.time()
    #
    if debug:
        print("Skipping full LAPACK diagonalization.")
        t2 = t1
    else:
        control(mat)
        t2 = time.time()
    #
    experiment(mat,vec,thresh,psize)
    t3 = time.time()
    #
    tmat,tdiag,tlanc,ttot = (t1-t0)/60 , (t2-t1)/60 , (t3-t2)/60 , (t3-t0)/60
    print("matrix build:          {} min".format(tmat))
    print("full diagonalization:  {} min".format(tdiag))
    print("lanczos extraction:    {} min".format(tlanc))
    print("Total:                 {} min".format(ttot))

if __name__ ==  "__main__":  main(sys.argv[1:])
