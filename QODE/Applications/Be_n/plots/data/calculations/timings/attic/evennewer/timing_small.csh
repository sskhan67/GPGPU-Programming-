#!/bin/csh
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

set counts = (  02 03 04 05 06 07 08 09 10 11 12 13 14  )

foreach c ( $counts )
  python3 run.py   $c 4.5   /morescratch/adutoi/programs/Qode/Applications/Be_n/dimer_H/6-31g/run/molecule_CI   16-115-550   load=u_4.5_16-115-550_t_1e-6.pickle
end
