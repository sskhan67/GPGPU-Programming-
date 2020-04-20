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

import sys
from   cProfile    import run
from   pstats      import Stats as stats
from   diag_random import main

debug = False
if not debug:
    sys.stdout = open(sys.argv[1]+'/stdout','w')
    sys.stderr = open(sys.argv[1]+'/stderr','w')
timingfilein  =   sys.argv[1]+'/timing.pstats'
timingfileout =   sys.argv[1]+'/timing.txt'

run('main(sys.argv[1:])',timingfilein)

f = open(timingfileout,'w')
timings = stats(timingfilein,stream=f)
timings.strip_dirs().sort_stats('cumulative','time').print_stats()
f.close()
