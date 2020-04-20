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
import numpy
from qode.math import *
from qode.atoms.integrals.basis import contractions, basis_function
from qode.atoms.integrals.integrals_better import kinetic, charge, charges, e_repulsion



N = charge(2,(0,0,0))
T = kinetic()
V = e_repulsion()

psi = basis_function(contractions["STO-3G"]["He"],(0,0,0))
psi /= sqrt(psi|psi)

E = 2*(psi|N|psi) + 2*(psi|T|psi) + ((psi,psi)|V|(psi,psi))
print(E)
