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
import qode.SCF_Code.main as scf_routine
import simple_cisd
import cisd_Lanczos

def print_amp( cisd_vec ):
    print("=======================================")
    print("-------------\nRef")
    print(cisd_vec.get_ref_amplitude())
    print("-------------\nSingles")
    print(cisd_vec.get_single_amplitude())
    print("-------------\nDoubles")
    print(cisd_vec.get_double_amplitude())
    print("=======================================")

#
# Energy calculated from Q-Chem 4.3
#
# SCF energy in the final basis set = -2.8551604262
# CCSD total energy                 = -2.87178917
# PSI4 CISD energy                  = -2.8717891505016
#

E_0, F_mat, V_mat, num_alpha_elec, num_beta_elecm, num_spatial_orb = scf_routine.scf_main("he.in", 1e-8)

print(E_0)




cisd_amp_vec = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elecm, num_spatial_orb )
H_cisd       = cisd_Lanczos.build_matrix( E_0, F_mat, V_mat )

evalue, evec = cisd_Lanczos.experiment(H_cisd, cisd_amp_vec, 1e-6)

print_amp(evec.vector)
print("EIGENVALUE =", evalue)
