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
import cisd_amplitude


class vector(space.vector):
    def __init__(self,vec):  
        self.vector = copy.deepcopy(vec)
        self.field = type(vec.get_coeffs()[0])

    @staticmethod
    def vecs_are_compatible(v,w):  
        return True

    @staticmethod
    def dot(v,w):
        # return a sum_i C_v(i) * C_w(i)
        # 
        # reference: 
        # print("type v =",type(v))
        # print("type w =",type(w))
        dot_prod = 1.0
        return dot_prod

    @staticmethod
    def sum(v,w):
        new_vec = vector(v.vector)
        return new_vec

    @staticmethod
    def scaled(c,v):
        new_vec = vector(v.vector)
        return new_vec  






class operator(space.operator):
    def __init__(self,mat): 
        self.H = mat
        self.field = type(self.H.get_hf_energy())

    @staticmethod
    def ops_are_compatible(A,B):  
        return True

    @staticmethod
    def op_and_vec_are_compatible(A,v):  
        return True

    @staticmethod
    def act_on_vec(A,v):
        return vector( A.H(v.vector) )

    @staticmethod
    def back_act_on_vec(v,A):
        # A is real symmetric!
        return self.act_on_vec(A,v)












