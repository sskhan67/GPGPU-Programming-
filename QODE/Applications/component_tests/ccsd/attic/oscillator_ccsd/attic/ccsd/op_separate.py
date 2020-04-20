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
from Applications.ccsd import noifccsd



def H1():
	pass

def H2():
	pass

def T1():
	pass

def T2():
	pass



def check_excitation_range(list_of_ops):
	#
	# Check H1 or H2
	#
	if H1 in list_of_ops:
		p = 1
	else:
		p = 2
	#
	# Get n_T
	#
	n_T = 0
	for item in list_of_ops:
		if item == T1:
			n_T += 1
		elif item == T2:
			n_T += 2
	#
	# Get k
	#
	k = len(list_of_ops) - 1
	#
	#
	low_limit  = n_T - p
	high_limit = n_T + p - k
	#
	#
	#
	if low_limit <= high_limit:
		if low_limit > 2 or high_limit < 1:
			return False
		else:
			# print(list_of_ops,"Limits = [%d, %d] n_T = %d  p = %d k = %d" %(low_limit,high_limit,n_T,p,k) )
			return True
	else:
		return False



def get_non_zero_operators(H1, H2, T1, T2):
	op_h = noifccsd.operator_list([[H1],[H2]], [1.0,1.0]) 
	op_t = noifccsd.operator_list([[T1],[T2]], [1.0,1.0]) 
	# print("H = H1 + H2")
	# op_h.print_info()
	# print("T = T1 + T2")
	# op_t.print_info()
	# term_count = 0
	op1  = noifccsd.commute(op_h, op_t)
	# term_count += len(op1.get_coeffs())
	all_op_terms = noifccsd.operator_list( copy.deepcopy(op1.get_operators()), copy.deepcopy(op1.get_coeffs()) )
	# print("Commute Once")
	# op1.print_info()
	for i in range(3):
		# print("Commute %d times" %(i+2))
		op1 = noifccsd.commute(op1, op_t)
		# term_count += len(op1.get_coeffs())
		all_op_terms.update_operators( all_op_terms.get_operators() + op1.get_operators() )
		new_coeffs = copy.deepcopy(op1.get_coeffs())
		if i == 0:
			new_coeffs = noifccsd.scale_py_list(new_coeffs, 0.5)
		elif i == 1:
			new_coeffs = noifccsd.scale_py_list(new_coeffs, 1.0/6.0)
		elif i == 2:
			new_coeffs = noifccsd.scale_py_list(new_coeffs, 1.0/24.0)
		all_op_terms.update_coeffs( all_op_terms.get_coeffs() + new_coeffs )
		# op1.print_info()

	
	operators = all_op_terms.get_operators()
	coeffs    = all_op_terms.get_coeffs()
	non_zero_ops = []
	op_coeffs = []
	for i in range(len(operators)):
		if check_excitation_range(operators[i]):
			non_zero_ops += [copy.deepcopy(operators[i])]
			op_coeffs    += [coeffs[i]]

	print("TOTAL NUM TERMS =", len(all_op_terms.get_operators()))
	print("NON-ZERO  TERMS =", len(non_zero_ops))
	return noifccsd.operator_list(non_zero_ops, op_coeffs)

if __name__ == "__main__":
	ops_obj = get_non_zero_operators(H1,H2,T1,T2)
	ops_obj.print_info()










