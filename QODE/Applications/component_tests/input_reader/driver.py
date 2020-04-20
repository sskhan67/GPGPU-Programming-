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
#!/bin/env python3.3

from sys       import argv
from math      import sqrt
from qode.util import read_input

print()
print(open(argv[1]).read())
print("============================\n")
data = read_input.from_file(argv[1],namespace=globals())
if len(argv)>2:  data = read_input.from_argv(argv[2:],write_to=data)
print(data.number)
print(data.letter)
print(data.fnexp)
print()
