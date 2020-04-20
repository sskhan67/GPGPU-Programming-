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
import qode
import numpy
from debug_CI import H_param
from my_space import space_traits
from class_CI import mat_dot_vec

def experiment(mat,vec,thresh,psize):
    S = qode.math.linear_inner_product_space(space_traits())
    Svec = S.member(vec)
    Smat = S.lin_op(mat)
    evalue,evec = qode.math.lanczos.lowest_eigen(Smat,Svec,thresh,psize)
    print("Lowest eigenvalue from Lanczos extraction:    {}".format(evalue))


#instantiate your mat_vec object

n_states = H_param["states_per_oscillator"] ** len(H_param["mass"])
M = mat_dot_vec ( H_param )                                

#instantiate a "guess" vector of the form expected by the __call__ of a mat_vec object

I      = numpy.matrix(numpy.zeros( (n_states, 1) ))
V      = I[:,0]
V[0,0] = 1

experiment(M,V,1e-3,9)
