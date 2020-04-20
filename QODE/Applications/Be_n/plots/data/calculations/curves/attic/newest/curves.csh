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



foreach d ( 3.8  3.9 )
  python3 run.py   2 $d   /morescratch/adutoi/programs/Qode/Applications/Be_n/dimer_H/6-31g/run/molecule_CI   16-115-550   load=u_4.5_16-115-550_t_1e-6.pickle
end

foreach i ( 4 5 6 7 8 9 )
  foreach j ( 0 1 2 3 4 5 6 7 8 9 )
    set d = $i.$j
    python3 run.py   2 $d   /morescratch/adutoi/programs/Qode/Applications/Be_n/dimer_H/6-31g/run/molecule_CI   16-115-550   load=u_4.5_16-115-550_t_1e-6.pickle
  end
end



foreach d ( 3.8  3.9 )
  python3 run.py   3 $d   /morescratch/adutoi/programs/Qode/Applications/Be_n/dimer_H/6-31g/run/molecule_CI   16-115-550   load=u_4.5_16-115-550_t_1e-6.pickle
end

foreach i ( 4 5 6 7 8 9 )
  foreach j ( 0 1 2 3 4 5 6 7 8 9 )
    set d = $i.$j
    python3 run.py   3 $d   /morescratch/adutoi/programs/Qode/Applications/Be_n/dimer_H/6-31g/run/molecule_CI   16-115-550   load=u_4.5_16-115-550_t_1e-6.pickle
  end
end
