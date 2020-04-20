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

# Dimer X2-CCSD(1e-6)
foreach i ( 3 4 5 6 7 8 9 )
  foreach j ( 0 1 2 3 4 5 6 7 8 9 )
    set d = $i.$j
    python3 run_ccsd.py   2 $d   ../../../../dimer_H/6-31g/H/16-115-550/load=H:16-115-550:thresh=1e-6:4.5:u.pickle
  end
end

# Dimer X2-CCSD(1e-6,noCT)
foreach i ( 3 4 5 6 7 8 9 )
  foreach j ( 0 1 2 3 4 5 6 7 8 9 )
    set d = $i.$j
    python3 run_ccsd.py   2 $d   ../../../../dimer_H/6-31g/H/0-115-0/load=H:0-115-0:thresh=1e-6:4.5:u.pickle
  end
end

# Trimer X2-CCSD(1e-6)
foreach i ( 3 4 5 6 7 8 9 )
  foreach j ( 0 1 2 3 4 5 6 7 8 9 )
    set d = $i.$j
    python3 run_ccsd.py   3 $d   ../../../../dimer_H/6-31g/H/16-115-550/load=H:16-115-550:thresh=1e-6:4.5:u.pickle
  end
end

# Trimer X2-CCSDT(1e-6)
foreach i ( 3 4 5 6 7 8 9 )
  foreach j ( 0 1 2 3 4 5 6 7 8 9 )
    set d = $i.$j
    python3 run_fci.py   3 $d   ../../../../dimer_H/6-31g/H/16-115-550/load=H:16-115-550:thresh=1e-6:4.5:u.pickle
  end
end

# Trimer X2-CCSDT(1e-7)
foreach i ( 3 4 5 6 7 8 9 )
  foreach j ( 0 1 2 3 4 5 6 7 8 9 )
    set d = $i.$j
    python3 run_fci.py   3 $d   ../../../../dimer_H/6-31g/H/16-115-550/load=H:16-115-550:thresh=1e-7:4.5:u.pickle
  end
end
