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
#numpy.set_printoptions(threshold=sys.maxsize)
numpy.set_printoptions(threshold=numpy.nan,linewidth=1000)
import copy
from numpy import linalg



def gen_config(x, H_param):
	list_of_states = []
	while x > 0:
		list_of_states.append( x % H_param["states_per_oscillator"] )
		x = x // H_param["states_per_oscillator"]
	pad = len(H_param["mass"]) - len(list_of_states)
	list_of_states = list_of_states + [0]*pad
	return (list_of_states)




def num_of_config(y, states_per_oscillator):
	N = 0
	for i in range( len(y) ):
		n = y[i] * (states_per_oscillator ** i)
		N += n
	return(N)



def compare(bra,ket, mass):
	diff = []
	for i in range( len(mass) ):
		if bra[i] != ket[i]:
			diff.append(( i, bra[i], ket[i] ))
	return(diff)
 



def CI_prod_intern(bra,ket, H_param):
  n_oscillators = len( H_param["mass"] )
  mat_elem = 0
  diffs = compare( bra, ket, H_param["mass"])
  if len(diffs) ==  0:          #Evaluating the diagonal elements
   for k in range( n_oscillators ):
    ei = (ket[k] + 0.5) * H_param["hbar"] * numpy.sqrt(H_param["k_OHconstant"][k]/H_param["arbitrary_mass"] )
    mat_elem += ei
   return mat_elem
  if len(diffs) ==  2:          #Evaluating non-diagonal elements
   pos1 = diffs[0][0]
   pos2 = diffs[1][0]
   bra1 = diffs[0][1] 
   ket1 = diffs[0][2]
   bra2 = diffs[1][1]
   ket2 = diffs[1][2]
   if (pos2 - pos1 == 1) or (pos2 - pos1 == n_oscillators - 1 ):
    p = 0.0
    for m in range(2):
     if abs(bra1 - ket1) == 1 and abs(bra2 - ket2) == 1 :
      pos_m = diffs[m][0]
      bra_m = diffs[m][1]
      ket_m = diffs[m][2]
      if bra_m - ket_m == 1:
       pi = ( (H_param["hbar"] / (2 * (H_param["mass"][pos_m] * H_param["k_OHconstant"][pos_m]) ** 0.5)  ) ** 0.5) * ((ket_m + 1 )**0.5)
       kij_OH = pos1    
      if (bra_m - ket_m) == -1:
       pi = ( (H_param["hbar"]/ (2 * (H_param["mass"][pos_m] * H_param["k_OHconstant"][pos_m]) ** 0.5)  ) ** 0.5) * ((ket_m)**0.5)
       kij_OH = pos1    
      if mat_elem == 0:
       mat_elem = H_param["k_coupling"][kij_OH] * pi
      else: 
       mat_elem *= pi
  return mat_elem



def evaluate(bra_string, final_vec, i, v, ket_string, CI_mat_element, H_param):
	k             =  num_of_config( bra_string, H_param["states_per_oscillator"] )
	mat_elem      =  CI_mat_element( bra_string, ket_string)
	final_vec[i] +=  mat_elem * v[k,0]



def eval_up_up( up_lim, final_vec, i, vec, ket_string, j, CI_mat_element, H_param):
	bra_string = copy.deepcopy( ket_string )
	if bra_string[j] < up_lim    and   bra_string[j+1] < up_lim:
		bra_string[j]   = bra_string[j]    + 1
		bra_string[j+1] = bra_string[j+1]  + 1
		evaluate( bra_string, final_vec, i, vec, ket_string, CI_mat_element, H_param)




def eval_dwn_dwn( up_lim, final_vec, i, vec, ket_string, j, CI_mat_element, H_param):
	bra_string  =  copy.deepcopy( ket_string )
	if bra_string[j] > 0 and bra_string[j+1] > 0:
		bra_string[j]   =  bra_string[j]   - 1
		bra_string[j+1] =  bra_string[j+1] - 1
		evaluate(bra_string, final_vec, i, vec, ket_string, CI_mat_element, H_param)




def eval_dwn_up( up_lim, final_vec, i, vec, ket_string, j, CI_mat_element, H_param):
	bra_string  =  copy.deepcopy( ket_string )
	if bra_string[j] > 0   and   bra_string[j+1] < up_lim:
		bra_string[j]   =  bra_string[j]    - 1
		bra_string[j+1] =  bra_string[j+1]  + 1
		evaluate(bra_string, final_vec, i, vec, ket_string, CI_mat_element, H_param)




def eval_up_dwn( up_lim, final_vec, i, vec, ket_string, j, CI_mat_element, H_param):
	bra_string  =  copy.deepcopy( ket_string )
	if bra_string[j] < up_lim    and    bra_string[j+1] > 0:
		bra_string[j] = bra_string[j] + 1
		bra_string[j+1] = bra_string[j+1] - 1
		evaluate(bra_string, final_vec, i, vec, ket_string, CI_mat_element, H_param)



def eval_diagonal( up_lim, final_vec, i, vec, ket_string, CI_mat_element, H_param):
	bra_string = copy.deepcopy(ket_string)
	evaluate(bra_string, final_vec, i, vec, ket_string, CI_mat_element, H_param)


