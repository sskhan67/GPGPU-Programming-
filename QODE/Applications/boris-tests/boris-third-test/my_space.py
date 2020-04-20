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
from class_CI import mat_dot_vec
import numpy
from debug_CI import H_param
import copy

class space_traits(object):
    def __init__(self):
        self.field = numpy.float64

    @staticmethod
    def dot(v,w):
        v = numpy.array(v.T)[0]
        w = numpy.array(w.T)[0]
        return numpy.dot(w,v)      #<a float>

    @staticmethod
    def add_to(v,w,c=1):
        if c==1: v +=   w          # <add c*w to v in place>
        else:    v += c*w
    @staticmethod
    def scale(n,v):
        v *= n                     #<scale v by n in place>

    @staticmethod
    def copy(v):
        return copy.deepcopy(v)    #<return an independent copy of v ... just use deepcopy?>

    @staticmethod
    def act_on_vec(op,v):          #<implement me>  
        return op(v)


    # Don't worry about this stuff ...
    def check_member(self,v):
        pass
    def check_lin_op(self,op):
        return False
    @staticmethod
    def back_act_on_vec(v,op):
        raise NotImplementedError
    @staticmethod
    def diagonal(op):
        raise NotImplementedError
    @staticmethod
    def function_on_diags(func,op):
        raise NotImplementedError
