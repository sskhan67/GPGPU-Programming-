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
the minimum time for two trials of each calculation was taken as the timing

  1270	10:04	columns 100 trial-1/CCSD_timing.txt trial-2/CCSD_timing.txt | awk '{if($4<$12){print($4)}else{print($12)}}' > CCSD_timing.txt
  1273	10:05	columns 100 trial-1/_T_timing.txt   trial-2/_T_timing.txt   | awk '{if($4<$12){print($4)}else{print($12)}}' >   _T_timing.txt
