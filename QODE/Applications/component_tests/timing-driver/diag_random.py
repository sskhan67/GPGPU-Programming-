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
import qode

def main(argv):
    data = qode.util.read_input.from_file(argv[0]+'/params.py')
    mat = qode.util.random_matrix(data.dim)
    eigenvals,eigenvecs = numpy.linalg.eigh(mat)
    lowest = qode.util.min.index(eigenvals)
    evalue = eigenvals[lowest]
    print("Lowest eigenvalue from full diagonalization:  {}".format(evalue))

if __name__ ==  "__main__":  main(sys.argv[1:])
