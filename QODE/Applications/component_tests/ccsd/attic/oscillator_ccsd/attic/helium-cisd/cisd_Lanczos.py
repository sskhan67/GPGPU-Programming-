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
import simple_cisd
import cisd_amplitude
import normal_order_hamiltonian

# This is how you call the Lanczos extraction routine.
# 
# 
# def experiment(mat,vec,thresh):
#     H     = simple_cisd.operator(mat)
#     guess = simple_cisd.vector(vec)
#     dim = vec.shape[0]//2     # arbitrary
#     if dim>50:  
#         dim = 50  # arbitrary
#     else: 
#         dim=4
#     evalue,evec = qode.math.lanczos.lowest_eigen(H,guess,thresh,dim=dim)
#     return evalue
# 
# 
#  MOST IMPORTANTLY!
#  evalue,evec = qode.math.lanczos.lowest_eigen(H,guess,thresh,dim=dim)
# 
#



def build_vector( num_alpha_elec, num_beta_elecm, num_spatial_orb ):
    # This is going to build a new reference amplitude vector for CISD.
    #
    return cisd_amplitude.cisd_amplitude(num_alpha_elec, num_beta_elecm, num_spatial_orb)

def build_matrix( E_0, F_mat, V_mat ):
    # HF energy, Fock Matrix in spin orbitals, Repulsion Matrix in spin orbitals
    #
    return normal_order_hamiltonian.normal_order_hamiltonian( E_0, F_mat, V_mat )



def experiment(mat,vec,thresh):
    print("===================================================")
    print("Running Lanczos Extraction Routine")
    H     = simple_cisd.operator(mat)
    guess = simple_cisd.vector(vec)
    print("vector dimension =",vec.get_vec_dimension())
    # dim = vec.get_vec_dimension()//2     # arbitrary
    # if dim>50:  
    #     dim = 50  # arbitrary
    # else: 
    #     dim=4     # arbitrary number
    dim = 10
    print("Extraction Space Dimension =",dim)
    evalue,evec = qode.math.lanczos.lowest_eigen(H,guess,thresh,dim=dim)
    print("Lowest eigenvalue from Lanczos extraction:    {}".format(evalue))
    print("===================================================")
    return evalue, evec







if __name__ ==  "__main__":  
    main(sys.argv[1:])
