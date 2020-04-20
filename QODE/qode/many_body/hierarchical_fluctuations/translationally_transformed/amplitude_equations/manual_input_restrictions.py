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
from index_restriction import get_restrictions

# Input code

order           = list(input("Type the order of the indices: "))
total_operators = int(input("What is the total number of operators?: "))
all_operators   = []
for i in range(total_operators):
    operator_indices = "Type the indices of operator {}: ".format(i+1)
    operator_i       = list(input(operator_indices))
    all_operators.append(operator_i)

# Summations dependencies ( the output will be written in latex )

restrictions = get_restrictions(order, all_operators) 
for index in restrictions:
	print(index[0])   
