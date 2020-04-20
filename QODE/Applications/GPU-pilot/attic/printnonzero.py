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


class bigger(object):
	def __init__(self, thresh):  self.thresh = thresh
	def __call__(self, element): return (abs(element)>self.thresh)

def space(indices, digits):
	string = ""
	sep = "["
	for index in indices:
		string += "{}{{:{}d}}".format(sep, digits).format(index)
		sep = ","
	return string+"]"

def printnonzero(element, index_digits=3, element_format="{}", thresh_check=bigger(1e-6), _superindices=[]):
	try:
		for i,subelement in enumerate(element):  printnonzero(subelement, index_digits, element_format, thresh_check, _superindices+[i])
	except:
		if thresh_check(element):  print(("{}: "+element_format).format(space(_superindices,index_digits), element))



#testarray = [[1,1,1],[2,4,8],[3,9,27]]

#printnonzero(testarray, index_digits=1, element_format="{:5d}", thresh_check=bigger(2))
