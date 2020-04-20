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
Working Code:
	First Working CCSD Code:   commutator.py
	No if conditions Code:     noifccsd.py
	CCSD with operator arithmetics: ccsd_separate_ops.py
	CCSD with Quasi-Newton PT: ccsd_qn_pt_simple_H.py
	CCSD with Quasi-Newton PT and DIIS acceleration: ccsd_qn_diss_pt_simple_H.py 

h_mat.npy:  Kinetic + Nuclear Attraction Energy Matrix from uncontracted HF/6-31g calculation
V_mat.npy:  Repulsion Matrix from the above calculation. It's this <i|<j|V|k>|l>, NOT the <ij||kl>

Three scans will be carried out.

CCSD QN PT 
CCSD QN PT DIIS starting at Cycle One
CCSD QN PT DIIS starting at Cycle Five
