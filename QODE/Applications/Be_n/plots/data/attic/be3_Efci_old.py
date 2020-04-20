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
from data.extract import X, cast
from data import atom

def parse(n, fields):  return float(fields[0]), float(fields[-1])-3*atom.energy

data = X(
parse,
"""\
4.0 -43.8381230689
4.1 -43.8383763065
4.2 -43.8385601873
4.3 -43.838688038
4.4 -43.8387712108
4.5 -43.838819907
4.6 -43.8388405881
4.7 -43.8388405508
4.8 -43.8388251195
4.9 -43.8387986223
5.0 -43.8387645556
5.1 -43.8387257205
5.2 -43.8386843388
5.3 -43.8386421467
5.4 -43.8386003561
5.5 -43.8385602145
5.6 -43.8385222756
5.7 -43.8384870042
5.8 -43.8384546713
5.9 -43.8384253918
6.0 -43.8383991593
""")
