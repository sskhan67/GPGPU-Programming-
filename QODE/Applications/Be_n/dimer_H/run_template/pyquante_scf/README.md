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

1) To build shared libraries needed, run:
   ./install.sh

2) To generate and run all inputs, use:
   ./runall.csh <basis>
   Alternatively, for a single internuclear distance, run, forexample:
   ./scf.sh Be.in <basis>
   ./scf.sh Be2_5_0ang.in <basis>
   <basis> is 6-31g or aug-cc-pvtz.  All the matrices are saved after the calculation is done.

3) Since the atom files never change:
   mv Be_*_spin.npy ../molecule_CI

4) Because these depend on geometry:
   mv Be2_*_spin.npy ..

5) We are done with this folder. (Optionally you can run ./clean.sh to remove everything.)

