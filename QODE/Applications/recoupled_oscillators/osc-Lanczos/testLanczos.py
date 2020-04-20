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
from time import time
import numpy as np
import qode
import simple_fast
from qode.coupled_oscillators import nbody_hamiltonian_func
import simple
#import parser_two_group as parser
import parser_six_identical_inline as parser
from qode.coupled_oscillators.analytical import analytical_solution

def build_vector(hamil_dim):
    vector = np.matrix( np.zeros((hamil_dim,1)) )
    vector[0,0] = 1.0
    return vector 

def build_matrix(rec_calc_obj):
    return rec_calc_obj

def parse_input(filename):
    # Read input and execute it.
    data = qode.util.read_input.from_file(filename)
    
    # Pass the data to any specific parser aliased as "parser" (module)
    # Get two objects, one is molecule obj, the other is recoupled obj.
    mol_obj, rec_obj = parser.build_recouple_object( data )
    
    
    dim = rec_obj.get_hamiltonian_dimension()

    mat = build_matrix(rec_obj)
    vec = build_vector(dim)
    
    thresh = 1e-3
    #print("NO PROBLEM FROM THE END OF parse_input")
    return mat, vec, thresh, mol_obj

def control(rec_obj):  # This is Lapack Diagonalization ( calling through Numpy )
    
    full_mat = nbody_hamiltonian_func.rec_obj_build_full_hamiltonian( rec_obj )
    eigenvals,eigenvecs = np.linalg.eigh(full_mat)
    lowest = qode.util.min.index(eigenvals)
    evalue = eigenvals[lowest]
    U = np.matrix(eigenvecs)
    #
    print("Lowest eigenvalue from full diagonalization:  {}".format(evalue))
    #print(lowest)
    #diag = U.T * mat * U
    #print(diag)

def lanczos_on_full_mat(mat,vec,thresh):
    H     = simple.operator(mat)
    guess = simple.vector(vec)
    print("vector dimension =",vec.shape)
    dim = vec.shape[0]//2     # arbitrary
    print("space dimension =",dim)
    if dim>50:  dim = 50  # arbitrary
    else: dim=4
    print("NEW DIM =",dim)
    #print("Lanczos: using subspace of size {}.".format(dim))
    evalue,evec = qode.math.lanczos.lowest_eigen(H,guess,thresh,dim=dim)
    print("Lowest eigenvalue from full mat Lanczos extraction: {}".format(evalue))


def experiment(mat,vec,thresh):
    print("===================================================")
    print("From Lanczos Extraction Function")
    H     = simple_fast.operator(mat)
    guess = simple_fast.vector(vec)
    print("vector dimension =",vec.shape)
    dim = vec.shape[0]//2     # arbitrary
    print("space dimension =",dim)
    if dim>50:  dim = 50  # arbitrary
    else: dim=4
    print("NEW DIM =",dim)
    #print("Lanczos: using subspace of size {}.".format(dim))
    evalue,evec = qode.math.lanczos.lowest_eigen(H,guess,thresh,dim=dim)
    #print("Lowest eigenvalue from Lanczos extraction:    {}".format(evalue))
    print("===================================================")
    return evalue


def get_analytical_solution( mol_obj ):
    required_states = 5
    k_eff, states_configure, E_analytical, vectors = \
        analytical_solution.analytical_main(
                                            mol_obj.list_of_masses() ,
                                            mol_obj.get_k_mat()      ,
                                            required_states,
                                            1.0,
                                            1.0,
                                            1.0
                                            )
    #print("Lowest Analytical Solution = ", sorted( E_analytical )[0] )
    return sorted( E_analytical )[0]


def main(argv):
    # when have better IO system, pass file pointer rather than use print statements
    t0 = time()
    mat,vec,thresh,mol_obj = parse_input(argv[0]+'/params.py')
    t1 = time()
    #control(mat)
    t2 = time()
    e_lanczos = experiment(mat,vec,thresh)
    t3 = time()
    tmat,tdiag,tlanc = (t1-t0)/60 , (t2-t1)/60 , (t3-t2)/60

    t4 = time()
    e_ana = get_analytical_solution( mol_obj )
    t5 = time()
    tanaly = (t5-t4)/60.0
    print("==========================================================")
    print("RESULTS")
    print("Lanczos threshold =  %.1e" %(thresh) )
    print("Lowest eigenvalue from Lanczos    Extraction  = {}".format(e_lanczos))
    print("Lowest eigenvalue from Analytical Solution    = {}".format(e_ana) )
    print("Accuracy = ", abs(e_lanczos - e_ana)/e_ana )
    print("==========================================================")
    print("TIMING INFO")
    print("matrix build:                   {} min".format(tmat))
    #print("full diagonalization:           {} min".format(tdiag))
    print("lanczos extraction:             {} min".format(tlanc))
    print("analytical solution:            {} min".format(tanaly))
    print("==========================================================")



if __name__ ==  "__main__":  main(sys.argv[1:])
