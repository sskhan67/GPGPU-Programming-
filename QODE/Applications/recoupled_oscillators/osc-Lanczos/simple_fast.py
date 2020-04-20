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
import copy
import numpy as np
import qode
from qode.math import space
# from qode.coupled_oscillators import mat_vec_product
from qode.coupled_oscillators.c_mat_vec_prod_routine import c_mat_vec_product
#from qode.coupled_oscillators.c_mat_vec_prod_routine.omp_prod import c_mat_vec_product

class vector(space.vector):
    def __init__(self,vec):  
        self.vector = copy.deepcopy(vec)
        self.field = type(vec[0,0])
    @staticmethod
    def vecs_are_compatible(v,w):  
        return True
    @staticmethod
    def dot(v,w):
        return (v.vector.T * w.vector)[0,0]
    @staticmethod
    def sum(v,w):
        return vector(v.vector + w.vector)
    @staticmethod
    def scaled(c,v):
        return vector(c * v.vector)  

class operator(space.operator):
    def __init__(self,mat): 
        self.eigval_list  = mat.get_sorted_eigen_energy()
        self.coupling_mat_list = mat.get_pair_dipole_mat()
        self.field = type(self.eigval_list[0][0])

    @staticmethod
    def ops_are_compatible(A,B):  
        return True

    @staticmethod
    def op_and_vec_are_compatible(A,v):  
        return True

    @staticmethod
    def act_on_vec(A,v):
        # A is an operator class instance, v is a vector class instance ( assume COLUMN Vector )
        # Convert eigval_list to PYTHON list of Numpy Arrays.
        # Convert coupling_mat_list to PYTHON list of Numpy 2D Arrays.
        eigen_list       =  [ np.array(item) for item in A.eigval_list ]
        dipole_mat_list  =  [ np.array(item) for item in A.coupling_mat_list ]
        
        new_vec = c_mat_vec_product.mat_vec_product_main(
                                           eigen_list ,        \
                                           dipole_mat_list ,   \
                                           v.vector              )
        # old Python Routine
        #        new_vec = c_mat_vec_product.mat_vec_product_main(
        #                                    eigen_list ,        \
        #                                    dipole_mat_list ,   \
        #                                    v.vector              )
        
        new_vec = vector( new_vec )
        return new_vec

    @staticmethod
    def back_act_on_vec(v,A):
        # A is real symmetric!
        return self.act_on_vec(A,v)












