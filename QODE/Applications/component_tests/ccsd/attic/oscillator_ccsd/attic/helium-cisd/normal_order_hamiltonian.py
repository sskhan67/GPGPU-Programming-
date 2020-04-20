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
import numpy as np
import copy
import cisd_amplitude
 

'''Full Hamiltonian:
\hat{H} = E_0  +  \sum_{ia} F_{ia} i^+ a  + \sum_{ij} - F_{ij} j i^+  + \sum_{ab} F_{ab} a^+ b  +  \sum_{ai} F_{ai} a^+ i
        + \sum_{ijab} \frac{1}{4} V_{ijab} i^+ j^+ a b +   \sum_{ijab} \frac{1}{4} V_{abij}    a^+b^+ ij
        + \sum_{ijkl} \frac{1}{4} V_{ijkl} kl i^+ j^+  + \sum_{abcd} \frac{1}{4} V_{abcd} a^+ b^+ c d +  \sum_{ijab} \frac{1}{4} V_{ajbi} a^+ i j^+ b  
        + \sum_{ia} \bigg[ \sum_{jk} \frac{1}{2} V_{ijak} k j^+ + \sum_{bc} - \frac{1}{2} V_{ibac} b^+ c \bigg] i^+ a
        + \sum_{ai}  a^+ i  \bigg[ \sum_{jk} \frac{1}{2} V_{jaki} k j^+ + \sum_{bc} - \frac{1}{2} V_{baci} b^+ c \bigg]'''
    



# This module takes CISD amplitudes and act the normal order Hamiltonian onto it to generate the new CISD amplitudes.
# There are 14 pieces in this normal order Hamiltonian.

# Each piece generates a piece of update and return it back. These updates will be summed over in the main function.


CISD_REF = True
CISD_SINGLE = True
CISD_DOUBLE = True




def ham1( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    
    # 1) onto ref
    if CISD_REF:
        # if E_0 * cisd_amp_obj.get_ref_amplitude() != 0.0:
            # print("h1 Ref =", E_0 * cisd_amp_obj.get_ref_amplitude() )
        cisd_amp_update_obj.update_ref_amplitude( E_0 * cisd_amp_obj.get_ref_amplitude() )
    
    # 2) onto singles
    if CISD_SINGLE:
        cisd_amp_update_obj.update_single_amplitude( E_0 * cisd_amp_obj.get_single_amplitude() )
    
    # 3) ontot doubles
    if CISD_DOUBLE:
        cisd_amp_update_obj.update_double_amplitude( E_0 * cisd_amp_obj.get_double_amplitude() )
    
    return cisd_amp_update_obj




def ham2( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    
    # 1) onto ref ==> no effect

    # 2) onto singles

   
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb

    # build F_eff = F_me ( occ - virt Fock Matrix )
    # Upper Submatrix
    F_aa  = F_mat[:num_alpha_elec, num_alpha_elec:num_spatial_orb]
    F_ab  = F_mat[:num_alpha_elec, num_spatial_orb+num_beta_elec:2*num_spatial_orb]
    F_up = np.concatenate((F_aa,F_ab) , axis=1)
    # Lower Submatrix
    F_ba  = F_mat[num_spatial_orb:num_spatial_orb+num_beta_elec, num_alpha_elec:num_spatial_orb]
    F_bb  = F_mat[num_spatial_orb:num_spatial_orb+num_beta_elec, num_spatial_orb+num_beta_elec:2*num_spatial_orb]
    F_low = np.concatenate((F_ba,F_bb) , axis=1)
    # Combine Upper and Lower
    F_eff = np.concatenate((F_up,F_low), axis=0)
    

    new_ref_amp = 0.0
    single_amp_mat = cisd_amp_obj.get_single_amplitude()
    
    for m  in range( num_occ_orb ):
        for e in range( num_virt_orb ):
            new_ref_amp += single_amp_mat[m,e] * F_eff[m,e]
    
    if CISD_REF:
        # if new_ref_amp != 0.0:
            # print("h2 Ref =", new_ref_amp )
        cisd_amp_update_obj.update_ref_amplitude( new_ref_amp )
                

    # 3) ontot doubles
    # The equation can be rewritten into this form,
    # C'_m^e = - \sum_{n f, n>m, f>m} C_{mn}^{ef} [ F_ne + F_me + F_mf + F_nf ]
    #
    
    double_amp_mat = cisd_amp_obj.get_double_amplitude()
    new_single_amp_mat = np.zeros(cisd_amp_update_obj.get_single_amplitude().shape)
    
    #for i in range(new_double_amp_mat.shape[0]):
    #    new_double_amp_mat[i] = np.random.randint(0,100,new_double_amp_mat.shape[1])
    #print("NEW DOUBLE AMPLITUDES MATRIX")
    #print(new_double_amp_mat)
    #print(F_eff)
    # sub_mat_size = num_alpha_elec + num_beta_elec
    
    # m < n && e < f
    for m in range(num_occ_orb):
        for n in range(m+1, num_occ_orb):
            for e in range(num_virt_orb):
                for f in range(e+1, num_virt_orb):
                    common_factor = double_amp_mat[e*num_occ_orb + m, f*num_occ_orb + n]
                    new_single_amp_mat[m,f] +=  common_factor * F_eff[n,e]
                    new_single_amp_mat[n,f] += -common_factor * F_eff[m,e]
                    new_single_amp_mat[n,e] +=  common_factor * F_eff[m,f]
                    new_single_amp_mat[m,e] += -common_factor * F_eff[n,f]

    if CISD_SINGLE:                            
        cisd_amp_update_obj.update_single_amplitude( -1.0 *  new_single_amp_mat )
                            
    return cisd_amp_update_obj






def ham3( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref ==> No Effect
    # 2) onto singles
    
    # build F_eff = F_mj ( occ - occ Fock Matrix )

    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    # Upper Submatrix
    F_aa  = F_mat[:num_alpha_elec, :num_alpha_elec]
    F_ab  = F_mat[:num_alpha_elec, num_spatial_orb:num_spatial_orb + num_beta_elec]
    F_up = np.concatenate((F_aa,F_ab) , axis=1)
    # Lower Submatrix
    F_ba  = F_mat[num_spatial_orb:num_spatial_orb + num_beta_elec, :num_alpha_elec]
    F_bb  = F_mat[num_spatial_orb:num_spatial_orb + num_beta_elec, num_spatial_orb:num_spatial_orb + num_beta_elec  ]
    F_low = np.concatenate((F_ba,F_bb) , axis=1)
    # Combine Upper and Lower
    F_eff = np.concatenate((F_up,F_low), axis=0)
 
    
    single_amp_mat = cisd_amp_obj.get_single_amplitude()
    new_single_amp_mat = np.zeros( single_amp_mat.shape )
    

    for j in range(num_occ_orb):
        for e in range(num_virt_orb):
            for m in range(num_occ_orb):
                new_single_amp_mat[j,e]  -= single_amp_mat[m,e] * F_eff[m,j]
    
    #print(new_single_amp_mat)
    if CISD_SINGLE:
        cisd_amp_update_obj.update_single_amplitude( new_single_amp_mat )


                    
    # 3) ontot doubles
    #
                    
    double_amp_mat = cisd_amp_obj.get_double_amplitude()
    new_double_amp_mat = np.zeros( double_amp_mat.shape )
                    

    
    # sub_mat_size = num_alpha_elec + num_beta_elec
    
    for e in range( num_virt_orb ):
        for f in range( e+1, num_virt_orb ):
            for j in range(num_occ_orb):
                for m in range(j+1, num_occ_orb):
                    #print(e, f, j, m, e*sub_mat_size + j, f*sub_mat_size+m, new_double_amp_mat[e*sub_mat_size + j, f*sub_mat_size+m] )
                    for n in range(m+1, num_occ_orb):
                        new_double_amp_mat[e*num_occ_orb + j, f*num_occ_orb + m] += \
                             double_amp_mat[e*num_occ_orb + m, f*num_occ_orb + n] * F_eff[n,j]
                    for n in range(m):
                        new_double_amp_mat[e*num_occ_orb + j, f*num_occ_orb + m] += \
                            -double_amp_mat[e*num_occ_orb + n, f*num_occ_orb + m] * F_eff[n,j]
                    for n in range(j+1, num_occ_orb):
                        new_double_amp_mat[e*num_occ_orb + j, f*num_occ_orb + m] += \
                            -double_amp_mat[e*num_occ_orb + j, f*num_occ_orb + n] * F_eff[n,m]
                    for n in range(j):
                        new_double_amp_mat[e*num_occ_orb + j, f*num_occ_orb + m] += \
                             double_amp_mat[e*num_occ_orb + n, f*num_occ_orb + j] * F_eff[n,m]
     
    if CISD_DOUBLE:               
        cisd_amp_update_obj.update_double_amplitude( new_double_amp_mat )
                    

    return cisd_amp_update_obj





def ham4( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref ==> No Effect
    # 2) onto singles
    
    # build F_eff = F_ae ( virt - virt Fock Matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb

    # Upper Submatrix
    F_aa  = F_mat[num_alpha_elec:num_spatial_orb, num_alpha_elec:num_spatial_orb]
    F_ab  = F_mat[num_alpha_elec:num_spatial_orb, num_spatial_orb+num_beta_elec:2*num_spatial_orb]
    F_up = np.concatenate((F_aa,F_ab) , axis=1)
    # Lower Submatrix
    F_ba  = F_mat[num_spatial_orb+num_beta_elec:2*num_spatial_orb ,num_alpha_elec:num_spatial_orb]
    F_bb  = F_mat[num_spatial_orb+num_beta_elec:2*num_spatial_orb ,num_spatial_orb+num_beta_elec:2*num_spatial_orb]
    F_low = np.concatenate((F_ba,F_bb) , axis=1)
    # Combine Upper and Lower
    F_eff = np.concatenate((F_up,F_low), axis=0)
 
    
    single_amp_mat = cisd_amp_obj.get_single_amplitude()
    new_single_amp_mat = np.zeros( single_amp_mat.shape )
    
    for m in range(num_occ_orb):
        for a in range(num_virt_orb):
            for e in range(num_virt_orb):
                new_single_amp_mat[m,a] += single_amp_mat[m,e] * F_eff[a,e]
    
    if CISD_SINGLE:
        cisd_amp_update_obj.update_single_amplitude( new_single_amp_mat )
                    
    # 3) ontot doubles
    
    double_amp_mat = cisd_amp_obj.get_double_amplitude()
    new_double_amp_mat = np.zeros( double_amp_mat.shape )
    
    # sub_mat_size = num_alpha_elec + num_beta_elec
                    
    for a in range(num_virt_orb):
        for f in range(a+1, num_virt_orb):
            for m in range(num_occ_orb):
                for n in range(m+1, num_occ_orb):
                    for e in range(f):
                        new_double_amp_mat[a*num_occ_orb + m, f*num_occ_orb + n] += \
                             double_amp_mat[e*num_occ_orb + m, f*num_occ_orb + n] * F_eff[a,e]
                    for e in range(f+1,num_virt_orb):
                        new_double_amp_mat[a*num_occ_orb + m, f*num_occ_orb + n] += \
                            -double_amp_mat[f*num_occ_orb + m, e*num_occ_orb + n] * F_eff[a,e]
                    for e in range(a):
                        new_double_amp_mat[a*num_occ_orb + m, f*num_occ_orb + n] += \
                            -double_amp_mat[e*num_occ_orb + m, a*num_occ_orb + n] * F_eff[f,e]
                    for e in range(a+1,num_virt_orb):
                        new_double_amp_mat[a*num_occ_orb + m, f*num_occ_orb + n] += \
                             double_amp_mat[a*num_occ_orb + m, e*num_occ_orb + n] * F_eff[f,e]
    
    if CISD_DOUBLE:
        cisd_amp_update_obj.update_double_amplitude( new_double_amp_mat )
    
    return cisd_amp_update_obj




def ham5( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref
    
    # build F_eff = F_ae ( virt - occ Fock Matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    # Upper Submatrix
    F_aa  = F_mat[num_alpha_elec:num_spatial_orb, :num_alpha_elec]
    F_ab  = F_mat[num_alpha_elec:num_spatial_orb, num_spatial_orb:num_spatial_orb+num_beta_elec]
    F_up = np.concatenate((F_aa,F_ab) , axis=1)
    # Lower Submatrix
    F_ba  = F_mat[num_spatial_orb+num_beta_elec:2*num_spatial_orb, :num_alpha_elec]
    F_bb  = F_mat[num_spatial_orb+num_beta_elec:2*num_spatial_orb, num_spatial_orb:num_spatial_orb+num_beta_elec]
    F_low = np.concatenate((F_ba,F_bb) , axis=1)
    # Combine Upper and Lower
    F_eff = np.concatenate((F_up,F_low), axis=0)
    
    
    new_single_amp_mat = np.zeros( cisd_amp_update_obj.get_single_amplitude().shape )
    ref_amp = cisd_amp_obj.get_ref_amplitude()
    
    for i in range(num_occ_orb):
        for a in range(num_virt_orb):
            new_single_amp_mat[i,a] = ref_amp * F_eff[a,i]

    if CISD_SINGLE:
        cisd_amp_update_obj.update_single_amplitude( new_single_amp_mat )
    
    
                
    # 2) onto singles
    new_double_amp_mat = np.zeros( cisd_amp_obj.get_double_amplitude().shape )
    single_amp_mat = cisd_amp_obj.get_single_amplitude()
                
    # sub_mat_size = num_alpha_elec + num_beta_elec
                
    for a in range(num_virt_orb):
        for e in range(a+1, num_virt_orb):
            for i in range(num_occ_orb):
                for m in range(i+1, num_occ_orb):
                    new_double_amp_mat[a*num_occ_orb + i, e*num_occ_orb + m] = \
                        - single_amp_mat[m,e] * F_eff[a,i] + single_amp_mat[i,e] * F_eff[a,m] + \
                            single_amp_mat[m,a] * F_eff[e,i] - single_amp_mat[i,a] * F_eff[e,m] 

    if CISD_DOUBLE:
        cisd_amp_update_obj.update_double_amplitude( -1.0 *  new_double_amp_mat )
    
    
    # 3) ontot doubles ==> No Effect
    #
    # TRIPLES ARE NOT CONSIDERED IN THIS CISD SCOPE
    #

    return cisd_amp_update_obj




def ham6( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref  ==> No Effect
    
    # 2) onto singles  ==> No Effect
    
    # 3) ontot doubles
    
    # build V_eff = V_occ,occ,virt,virt ( V is the  repulsion matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_spin_orb       = 2 * num_spatial_orb
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    
    new_ref_amp = 0.0
    double_amp_mat = cisd_amp_obj.get_double_amplitude()
    
    
    index_i = [ i for i in range(num_alpha_elec)] + [ num_spatial_orb +  i for i in range(num_beta_elec) ]
    index_j = index_i
    index_k = [ num_alpha_elec + i for i in range(num_alpha_virt_orb)] + [ num_spatial_orb + num_beta_elec + i for i in range(num_beta_virt_orb) ]
    index_l = index_k
    
    
    J_eff = []
    
    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_k:
                for l in index_l:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            J_eff += [ V_line_tmp ]
                                    
    J_eff = np.matrix(J_eff)
    # print(J_eff)


                
    sub_mat_size = num_alpha_elec + num_beta_elec
    
    for e in range(num_virt_orb):
        for f in range(e+1, num_virt_orb):
            for m in range(num_occ_orb):
                for n in range(m+1, num_occ_orb):
                    new_ref_amp += double_amp_mat[e*sub_mat_size + m, f*sub_mat_size + n] * (-1.0) * \
                        ( J_eff[m*num_occ_orb + n, e*num_virt_orb + f ] - J_eff[m*num_occ_orb + n, f*num_virt_orb + e ] )
    # Four time of V_mnef cancel out the 1/4 in the front
    # print('h6 ref =', new_ref_amp)


    if CISD_REF:
        # if new_ref_amp != 0.0:
            # print("h6 Ref =", new_ref_amp)
        cisd_amp_update_obj.update_ref_amplitude( -1.0 * new_ref_amp )
                
    return cisd_amp_update_obj




def ham7( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref
    
    # build V_eff = V_virt,virt,occ,occ ( V is the  repulsion matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_spin_orb       = 2 * num_spatial_orb
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    ref_amp = cisd_amp_obj.get_ref_amplitude()  
    new_double_amp_mat = np.zeros( cisd_amp_obj.get_double_amplitude().shape )
  
    index_i = [ num_alpha_elec + i for i in range(num_alpha_virt_orb)] + [ num_spatial_orb + num_beta_elec + i for i in range(num_beta_virt_orb)]
    index_j = index_i
    index_k = [i for i in range(num_alpha_elec)] + [ num_spatial_orb +  i for i in range(num_beta_elec) ]
    index_l = index_k
  
    J_eff = []
    
    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_k:
                for l in index_l:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            J_eff += [ V_line_tmp ]
    
    J_eff = np.matrix(J_eff)
    
    sub_mat_size = num_alpha_elec + num_beta_elec
    
    for a in range(num_virt_orb):
        for b in range(a+1, num_virt_orb):
            for i in range(num_occ_orb):
                for j in range(i+1, num_occ_orb):
                    new_double_amp_mat[ a*sub_mat_size + i, b*sub_mat_size + j ] =  ref_amp * (-1.0) *\
                        ( J_eff[ a* num_virt_orb + b, i*num_occ_orb + j ] -  J_eff[ a* num_virt_orb + b, j*num_occ_orb + i ] )
    # Four time of V_mnef cancel out the 1/4 in the front
    
    if CISD_DOUBLE:
        cisd_amp_update_obj.update_double_amplitude( -1.0 * new_double_amp_mat )
    
    
    # 2) onto singles  ==> No Effect
    
    # 3) ontot doubles ==> No Effect
                            
                            
    return cisd_amp_update_obj




def ham8( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref  ==> No Effect
    
    # 2) onto singles  ==> No Effect
    
    # 3) ontot doubles
    
    # build V_eff = V_occ,occ,occ,occ ( V is the  repulsion matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_spin_orb       = 2 * num_spatial_orb
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb

    
    index_i = [ i for i in range(num_alpha_elec)] + [ num_spatial_orb + i for i in range(num_beta_elec) ]
    index_j = index_i
    index_k = index_i
    index_l = index_i
    
    J_eff = []
    
    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_k:
                for l in index_l:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            J_eff += [ V_line_tmp ]
    
    J_eff = np.matrix(J_eff)
    
    sub_mat_size = num_alpha_elec + num_beta_elec

    double_amp_mat = cisd_amp_obj.get_double_amplitude()
    new_double_amp_mat = np.zeros( double_amp_mat.shape )
    
    for e in range(num_virt_orb):
        for f in range(e+1, num_virt_orb):
            for k in range(num_occ_orb):
                for l in range(k+1, num_occ_orb):
                    for m in range(num_occ_orb):
                        for n in range(m+1, num_occ_orb):
                            new_double_amp_mat[ e*sub_mat_size + k, f*sub_mat_size + l ] +=  \
                                 double_amp_mat[ e*sub_mat_size + m, f*sub_mat_size + n ] * \
                                ( J_eff[ m*num_occ_orb + n, k*num_occ_orb + l ] - J_eff[ m*num_occ_orb + n, l*num_occ_orb + k ] )   
    # Four time of V_mnef cancel out the 1/4 in the front
    if CISD_DOUBLE:
        cisd_amp_update_obj.update_double_amplitude( new_double_amp_mat )
    
    
    return cisd_amp_update_obj




def ham9( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref  ==> No Effect
    
    # 2) onto singles  ==> No Effect
    
    # 3) ontot doubles
    
    # build V_eff = V_virt,virt,virt,virt ( V is the  repulsion matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_spin_orb       = 2 * num_spatial_orb
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    
    index_i = [ num_alpha_elec + i for i in range(num_alpha_virt_orb)] + [ num_spatial_orb + num_beta_elec + i for i in range(num_beta_virt_orb) ]
    index_j = index_i
    index_k = index_i
    index_l = index_i

    J_eff = []
    
    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_k:
                for l in index_l:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            J_eff += [ V_line_tmp ]
    
    J_eff = np.matrix(J_eff)
    
    
    sub_mat_size = num_alpha_elec + num_beta_elec
    
    double_amp_mat = cisd_amp_obj.get_double_amplitude()
    new_double_amp_mat = np.zeros( double_amp_mat.shape )
    
    
    for a in range(num_virt_orb):
        for b in range(a+1, num_virt_orb):
            for m in range(num_occ_orb):
                for n in range(m+1, num_occ_orb):
                    for e in range(num_virt_orb):
                        for f in range(e+1, num_virt_orb):
                            new_double_amp_mat[a*sub_mat_size+m, b*sub_mat_size+n] += \
                                double_amp_mat[e*sub_mat_size+m, f*sub_mat_size+n] * \
                                ( J_eff[a*num_virt_orb + b, e*num_virt_orb + f] - J_eff[a*num_virt_orb + b, f*num_virt_orb + e] )
    # FOUR TIMES the V_abef cancel out the 1/4
    
    if CISD_DOUBLE:
        cisd_amp_update_obj.update_double_amplitude( new_double_amp_mat )
    
    return cisd_amp_update_obj




def ham10( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref  ==> No Effect

    
    # 2) onto singles
    
    # build V_eff = V_virt,occ,virt,occ ( V is the repulsion matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_spin_orb       = 2 * num_spatial_orb
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    
    index_i = [ num_alpha_elec + i for i in range(num_alpha_virt_orb)] + [ num_spatial_orb + num_beta_elec + i for i in range(num_beta_virt_orb) ]
    index_j = [ i for i in range(num_alpha_elec) ] + [ num_spatial_orb + i for i in range(num_beta_elec) ]
    index_k = index_i
    index_l = index_j

    # Coulomb Matrix
    J_eff = []
    
    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_k:
                for l in index_l:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            J_eff += [ V_line_tmp ]
    
    # Exchange Matrix 
    K_eff = []

    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_l:
                for l in index_k:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            K_eff += [ V_line_tmp ]


    J_eff = np.matrix(J_eff)
    K_eff = np.matrix(K_eff)
    
    single_amp_mat = cisd_amp_obj.get_single_amplitude()
    new_single_amp_mat = np.zeros( single_amp_mat.shape )
                
    for i in range(num_occ_orb):
        for a in range(num_virt_orb):
            for m in range(num_occ_orb):
                for e in range(num_virt_orb):
                    new_single_amp_mat[i,a] += 0.25 * single_amp_mat[m,e] * (-1.0) * \
                        ( J_eff[a*num_occ_orb+m, e*num_occ_orb+i] - K_eff[a*num_occ_orb+m, i*num_virt_orb+e] )
    
    if CISD_SINGLE:
        cisd_amp_update_obj.update_single_amplitude( new_single_amp_mat )
    
    # 3) ontot doubles
    
    double_amp_mat = cisd_amp_obj.get_double_amplitude()
    new_double_amp_mat = np.zeros( double_amp_mat.shape )

    # C'_im^ae = C~_im^ae - C~_im^ea - C~_mi^ae + C~_mi^ea 
    
    # Loop start with,
    #   n < m ______ f < e
    #          |____ f > e
    #          |____ f < a
    #          |____ f > a
    #
    #   n > m ______ f < e
    #          |____ f > e
    #          |____ f < a
    #          |____ f > a
    #
    #   n < i ______ f < e
    #          |____ f > e
    #          |____ f < a
    #          |____ f > a
    #
    #   n > i ______ f < e
    #          |____ f > e
    #          |____ f < a
    #          |____ f > a
    #
    #
    sub_mat_size = num_alpha_elec + num_beta_elec
                            



    for a in range(num_virt_orb):
        for e in range(a+1, num_virt_orb):
            for i in range(num_occ_orb):
                for m in range(i+1, num_occ_orb):
                    # Start Inner Loop 
                    for n in range(m):
                        for f in range(e):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] +=\
                                    0.25* double_amp_mat[f*sub_mat_size+n, e*sub_mat_size+m] * (-1.0) * \
                                        ( J_eff[a*num_occ_orb+n, f*num_occ_orb+i] - K_eff[a*num_occ_orb+n, i*num_virt_orb+f] )
                        for f in range(e+1, num_virt_orb):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] -=\
                                    0.25* double_amp_mat[e*sub_mat_size+n, f*sub_mat_size+m] * (-1.0) * \
                                        ( J_eff[a*num_occ_orb+n, f*num_occ_orb+i] - K_eff[a*num_occ_orb+n, i*num_virt_orb+f] ) 
                        for f in range(a):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] -=\
                                    0.25* double_amp_mat[f*sub_mat_size+n, a*sub_mat_size+m] * (-1.0) * \
                                        ( J_eff[e*num_occ_orb+n, f*num_occ_orb+i] - K_eff[e*num_occ_orb+n, i*num_virt_orb+f] )
                        for f in range(a+1, num_virt_orb):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] +=\
                                    0.25* double_amp_mat[a*sub_mat_size+n, f*sub_mat_size+m] * (-1.0) * \
                                        ( J_eff[e*num_occ_orb+n, f*num_occ_orb+i] - K_eff[e*num_occ_orb+n, i*num_virt_orb+f] )

                    for n in range(m+1, num_occ_orb):
                        for f in range(e):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] -=\
                                    0.25* double_amp_mat[f*sub_mat_size+m, e*sub_mat_size+n] * (-1.0) * \
                                        ( J_eff[a*num_occ_orb+n, f*num_occ_orb+i] - K_eff[a*num_occ_orb+n, i*num_virt_orb+f] )
                        for f in range(e+1, num_virt_orb):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] +=\
                                    0.25* double_amp_mat[e*sub_mat_size+m, f*sub_mat_size+n] * (-1.0) * \
                                        ( J_eff[a*num_occ_orb+n, f*num_occ_orb+i] - K_eff[a*num_occ_orb+n, i*num_virt_orb+f] )
                        for f in range(a):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] +=\
                                    0.25* double_amp_mat[f*sub_mat_size+m, a*sub_mat_size+n] * (-1.0) * \
                                        ( J_eff[e*num_occ_orb+n, f*num_occ_orb+i] - K_eff[e*num_occ_orb+n, i*num_virt_orb+f] )
                        for f in range(a+1, num_virt_orb):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] -=\
                                    0.25* double_amp_mat[a*sub_mat_size+m, f*sub_mat_size+n] * (-1.0) * \
                                        ( J_eff[e*num_occ_orb+n, f*num_occ_orb+i] - K_eff[e*num_occ_orb+n, i*num_virt_orb+f] )

                    for n in range(i):
                        for f in range(e):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] -=\
                                    0.25* double_amp_mat[f*sub_mat_size+n, e*sub_mat_size+i] * (-1.0) * \
                                        ( J_eff[a*num_occ_orb+n, f*num_occ_orb+m] - K_eff[a*num_occ_orb+n, m*num_virt_orb+f] )
                        for f in range(e+1, num_virt_orb):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] +=\
                                    0.25* double_amp_mat[e*sub_mat_size+n, f*sub_mat_size+i] * (-1.0) * \
                                        ( J_eff[a*num_occ_orb+n, f*num_occ_orb+m] - K_eff[a*num_occ_orb+n, m*num_virt_orb+f] )
                        for f in range(a):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] +=\
                                    0.25* double_amp_mat[f*sub_mat_size+n, a*sub_mat_size+i] * (-1.0) * \
                                        ( J_eff[e*num_occ_orb+n, f*num_occ_orb+m] - K_eff[e*num_occ_orb+n, m*num_virt_orb+f] )
                        for f in range(a+1, num_virt_orb):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] -=\
                                    0.25* double_amp_mat[a*sub_mat_size+n, f*sub_mat_size+i] * (-1.0) * \
                                        ( J_eff[e*num_occ_orb+n, f*num_occ_orb+m] - K_eff[e*num_occ_orb+n, m*num_virt_orb+f] )

                    for n in range(i+1, num_occ_orb):
                        for f in range(e):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] +=\
                                    0.25* double_amp_mat[f*sub_mat_size+i, e*sub_mat_size+n] * (-1.0) * \
                                        ( J_eff[a*num_occ_orb+n, f*num_occ_orb+m] - K_eff[a*num_occ_orb+n, m*num_virt_orb+f] )
                        for f in range(e+1, num_virt_orb):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] -=\
                                    0.25* double_amp_mat[e*sub_mat_size+i, f*sub_mat_size+n] * (-1.0) * \
                                        ( J_eff[a*num_occ_orb+n, f*num_occ_orb+m] - K_eff[a*num_occ_orb+n, m*num_virt_orb+f] )
                        for f in range(a):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] -=\
                                    0.25* double_amp_mat[f*sub_mat_size+i, a*sub_mat_size+n] * (-1.0) * \
                                        ( J_eff[e*num_occ_orb+n, f*num_occ_orb+m] - K_eff[e*num_occ_orb+n, m*num_virt_orb+f] )
                        for f in range(a+1, num_virt_orb):
                            new_double_amp_mat[a*sub_mat_size+i, e*sub_mat_size+m] +=\
                                    0.25* double_amp_mat[a*sub_mat_size+i, f*sub_mat_size+n] * (-1.0) * \
                                        ( J_eff[e*num_occ_orb+n, f*num_occ_orb+m] - K_eff[e*num_occ_orb+n, m*num_virt_orb+f] )
                        
                        
    

    if CISD_DOUBLE:
        cisd_amp_update_obj.update_double_amplitude( new_double_amp_mat )

    return cisd_amp_update_obj



def ham11( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref  ==> No Effect
    # 2) onto singles  ==> No Effect
    
    # 3) ontot doubles
    
    # build V_eff = V_occ,occ,virt,occ ( V is the repulsion matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_spin_orb       = 2 * num_spatial_orb
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    
    index_i = [ i for i in range(num_alpha_elec) ] + [ num_spatial_orb + i for i in range(num_beta_elec) ]
    index_j = index_i
    index_k = [ num_alpha_elec + i for i in range(num_alpha_virt_orb)] + [ num_spatial_orb + num_beta_elec + i for i in range(num_beta_virt_orb) ]
    index_l = index_i
    
    J_eff = []
    
    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_k:
                for l in index_l:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            J_eff += [ V_line_tmp ]
    

    K_eff = []

    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_l:
                for l in index_k:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            K_eff += [ V_line_tmp ]


    J_eff = np.matrix(J_eff)
    K_eff = np.matrix(K_eff)    


    double_amp_mat = cisd_amp_obj.get_double_amplitude()
    new_single_amp_mat = np.zeros( cisd_amp_obj.get_single_amplitude().shape )
                
    sub_mat_size = num_alpha_elec + num_beta_elec

    for k in range(num_occ_orb):
        for e in range(num_virt_orb):
            for m in range(num_occ_orb):
                for n in range(m+1, num_occ_orb):
                    for f in range(e+1, num_virt_orb):
                        new_single_amp_mat[k,e]  += double_amp_mat[ e*sub_mat_size+m, f*sub_mat_size+n ] * (-1.0) *  \
                            ( J_eff[ m*num_occ_orb+n, f*num_occ_orb+k ] - K_eff[ m*num_occ_orb+n, k*num_virt_orb+f ] )
                    for f in range(e):
                        new_single_amp_mat[k,e]  += - double_amp_mat[ f*sub_mat_size+m, e*sub_mat_size+n ] * (-1.0) *  \
                            ( J_eff[ m*num_occ_orb+n, f*num_occ_orb+k ] - K_eff[ m*num_occ_orb+n, k*num_virt_orb+f ]  )
    
    if CISD_SINGLE:
        cisd_amp_update_obj.update_single_amplitude( new_single_amp_mat )

    return cisd_amp_update_obj




def ham12( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref  ==> No Effect
    # 2) onto singles  ==> No Effect
    
    # 3) ontot doubles
    
    # build V_eff = V_occ,virt,virt,virt ( V is the repulsion matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_spin_orb       = 2 * num_spatial_orb
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    
    index_i = [ i for i in range(num_alpha_elec) ] + [ num_spatial_orb + i for i in range(num_beta_elec) ]
    index_j = [ num_alpha_elec + i for i in range(num_alpha_virt_orb)] + [ num_spatial_orb + num_beta_elec + i for i in range(num_beta_virt_orb) ]
    index_k = index_j
    index_l = index_j
    
    J_eff = []
    
    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_k:
                for l in index_l:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            J_eff += [ V_line_tmp ]
    
    J_eff = np.matrix(J_eff)
    
    
    double_amp_mat = cisd_amp_obj.get_double_amplitude()
    new_single_amp_mat = np.zeros( cisd_amp_obj.get_single_amplitude().shape )

    sub_mat_size = num_alpha_elec + num_beta_elec

    for n in range(num_occ_orb):
        for b in range(num_virt_orb):
            for e in range(num_virt_orb):
                for f in range(e+1, num_virt_orb):
                    for m in range(num_occ_orb):
                        for n in range(m+1,num_occ_orb):
                            new_single_amp_mat[n,b] += double_amp_mat[ e*sub_mat_size+m, f*sub_mat_size+n ] * (-1.0) *  \
                                ( J_eff[ m*num_virt_orb+b, e*num_virt_orb+f ] - J_eff[ m*num_virt_orb+b, f*num_virt_orb+e ] )
                        for n in range(m):
                            new_single_amp_mat[n,b] += - double_amp_mat[ e*sub_mat_size+n, f*sub_mat_size+m ] * (-1.0) *  \
                                ( J_eff[ m*num_virt_orb+b, e*num_virt_orb+f ] - J_eff[ m*num_virt_orb+b, f*num_virt_orb+e ] )

    if CISD_SINGLE:
        cisd_amp_update_obj.update_single_amplitude( new_single_amp_mat )
    
    return cisd_amp_update_obj




def ham13( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    # 1) onto ref  ==> No Effect

    # 2) onto singles
    # build V_eff = V_occ,virt,occ,occ ( V is the repulsion matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_spin_orb       = 2 * num_spatial_orb
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    
    index_i = [ i for i in range(num_alpha_elec) ] + [ num_spatial_orb + i for i in range(num_beta_elec) ]
    index_j = [ num_alpha_elec + i for i in range(num_alpha_virt_orb)] + [ num_spatial_orb + num_beta_elec + i for i in range(num_beta_virt_orb) ]
    index_k = index_i
    index_l = index_i
    
    J_eff = []
    
    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_k:
                for l in index_l:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            J_eff += [ V_line_tmp ]
    
    J_eff = np.matrix(J_eff)


    single_amp_mat = cisd_amp_obj.get_single_amplitude()
    new_double_amp_mat = np.zeros( cisd_amp_obj.get_double_amplitude().shape )
    sub_mat_size = num_alpha_elec + num_beta_elec
    
    for a in range(num_virt_orb):
        for e in range(a+1, num_virt_orb):
            for i in range(num_occ_orb):
                for k in range(i+1, num_occ_orb):
                    for m in range(num_occ_orb):
                        new_double_amp_mat[ a*sub_mat_size+i, e*sub_mat_size+k ] +=  \
                                 single_amp_mat[m,e] * (-1.0) * \
                                 ( J_eff[m*num_virt_orb+a, i*num_occ_orb+k] - J_eff[m*num_virt_orb+a, k*num_occ_orb+i] ) \
                                 - single_amp_mat[m,a] * (-1.0) * \
                                 ( J_eff[m*num_virt_orb+e, i*num_occ_orb+k] - J_eff[m*num_virt_orb+e, k*num_occ_orb+i] )
    #  C'_ik^ae = \sum_m (C_me V_maik - C_ma V_meik)

    if CISD_DOUBLE:
        cisd_amp_update_obj.update_double_amplitude(  -1.0 * new_double_amp_mat )
    
    # 3) ontot doubles  ==> No Effect
    return cisd_amp_update_obj




def ham14( cisd_amp_obj , E_0 , F_mat , V_mat ):
    cisd_amp_update_obj = copy.deepcopy( cisd_amp_obj )
    cisd_amp_update_obj.clean_all_amplitude()
    
    # 1) onto ref  ==> No Effect
    
    # 2) onto singles
    
    # build V_eff = V_virt,virt,virt,occ ( V is the repulsion matrix )
    num_alpha_elec = cisd_amp_obj.get_num_alpha_elec()
    num_beta_elec  = cisd_amp_obj.get_num_beta_elec()
    num_alpha_virt_orb = cisd_amp_obj.get_num_alpha_virt_orb()
    num_beta_virt_orb  = cisd_amp_obj.get_num_bete_virt_orb()
    num_spatial_orb    = cisd_amp_obj.get_num_spatial_orb()
    num_spin_orb       = 2 * num_spatial_orb
    num_occ_orb  = num_alpha_elec + num_beta_elec
    num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb
    
    
    index_i = [ num_alpha_elec + i for i in range(num_alpha_virt_orb)] + [ num_spatial_orb + num_beta_elec + i for i in range(num_beta_virt_orb) ]
    index_j = index_i
    index_k = index_i
    index_l = [ i for i in range(num_alpha_elec) ] + [ num_spatial_orb + i for i in range(num_beta_elec) ]
    
    J_eff = []
    
    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_k:
                for l in index_l:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            J_eff += [ V_line_tmp ]
    
    J_eff = np.matrix(J_eff)
 

    K_eff = []

    for i in index_i:
        for j in index_j:
            V_line_tmp = []
            for k in index_l:
                for l in index_k:
                    V_line_tmp +=  [ V_mat[ i*num_spin_orb + j, k*num_spin_orb + l ] ]
            K_eff += [ V_line_tmp ]
    
    K_eff = np.matrix(K_eff)

    # print(J_eff)
    # print(K_eff)


    single_amp_mat = cisd_amp_obj.get_single_amplitude()
    new_double_amp_mat = np.zeros( cisd_amp_obj.get_double_amplitude().shape )
    sub_mat_size = num_alpha_elec + num_beta_elec

    for a in range(num_virt_orb):
        for b in range(a+1, num_virt_orb):
            for i in range(num_occ_orb):
                for m in range(i+1, num_occ_orb):
                    for e in range(num_virt_orb):
                        new_double_amp_mat[ a*sub_mat_size+i, b*sub_mat_size+m ] += single_amp_mat[i,e] * (-1.0) *  \
                        ( J_eff[a*num_virt_orb+b, e*num_occ_orb+m] - K_eff[a*num_virt_orb+b, m*num_virt_orb+e] ) \
                        - single_amp_mat[m,e] * (-1.0) * \
                        ( J_eff[a*num_virt_orb+b, e*num_occ_orb+i] - K_eff[a*num_virt_orb+b, i*num_virt_orb+e] )

    # C'_im^ab = \sum_e ( C_ie V_abem - C_me V_abei )

    if CISD_DOUBLE:
        cisd_amp_update_obj.update_double_amplitude( -1.0 * new_double_amp_mat )
                
    
    # 3) ontot doubles  ==> No Effect

    
    return cisd_amp_update_obj





def print_amp( cisd_vec ):
    print("=======================================")
    print("-------------\nRef")
    print(cisd_vec.get_ref_amplitude())
    print("-------------\nSingles")
    print(cisd_vec.get_single_amplitude())
    print("-------------\nDoubles")
    print(cisd_vec.get_double_amplitude())
    print("=======================================")


def add_and_update_cisd_amp_obj( destination_amp_obj, incremental_amp_obj  ):
    destination_amp_obj.update_ref_amplitude( destination_amp_obj.get_ref_amplitude() + incremental_amp_obj.get_ref_amplitude() )
    destination_amp_obj.update_single_amplitude( destination_amp_obj.get_single_amplitude() + incremental_amp_obj.get_single_amplitude() )
    destination_amp_obj.update_double_amplitude( destination_amp_obj.get_double_amplitude() + incremental_amp_obj.get_double_amplitude() )
    return destination_amp_obj



#===================== MAIN FUNCTION RESIDES HERE =========================================================



def normal_order_hamiltonian_act_on_vec( cisd_amp_obj, E_0, F_mat, V_mat ):
    # E_0 is the HF ground state energy
    # F_mat is the Fock matrix in spin orbital basis
    # V_mat is the Repulsion Matrix in spin orbital basis.
    cisd_new_amp_obj = copy.deepcopy( cisd_amp_obj )
    cisd_new_amp_obj.clean_all_amplitude()
    F_mat = np.matrix( F_mat )
    V_mat = np.matrix( V_mat )

    update_obj = ham1( cisd_amp_obj  , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 1")
    # print_amp(update_obj)
    
    update_obj = ham2( cisd_amp_obj  , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 2")
    # print_amp(update_obj)
            
    update_obj = ham3( cisd_amp_obj  , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 3")
    # print_amp(update_obj)

    update_obj = ham4( cisd_amp_obj  , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 4")
    # print_amp(update_obj)

    update_obj = ham5( cisd_amp_obj  , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 5")
    # print_amp(update_obj)
 
    update_obj = ham6( cisd_amp_obj  , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 6")
    # print_amp(update_obj)
    
    update_obj = ham7( cisd_amp_obj  , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 7")
    # print_amp(update_obj)

    update_obj = ham8( cisd_amp_obj  , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 8")
    # print_amp(update_obj)
    
    update_obj = ham9( cisd_amp_obj  , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 9")
    # print_amp(update_obj)
    
    update_obj = ham10( cisd_amp_obj , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 10")
    # print_amp(update_obj)   
    
    update_obj = ham11( cisd_amp_obj , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 11")
    # print_amp(update_obj)

    
    update_obj = ham12( cisd_amp_obj , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )
    
    # print("Part 12")
    # print_amp(update_obj)
    

    update_obj = ham13( cisd_amp_obj , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )

    # print("Part 13")
    # print_amp(update_obj)            

    update_obj = ham14( cisd_amp_obj , E_0 , F_mat , V_mat )
    add_and_update_cisd_amp_obj( cisd_new_amp_obj , update_obj )

    # print("Part 14")
    # print_amp(update_obj)
    # print("NEW AMP OBJ")
    # print_amp( cisd_new_amp_obj )
    # if cisd_new_amp_obj.get_ref_amplitude() != 0.0:
    #     print(cisd_new_amp_obj.get_ref_amplitude())

    return cisd_new_amp_obj





class normal_order_hamiltonian:
    def __init__(self, E_0, F_mat, V_mat):
        self.E_0   = np.matrix( np.zeros((1,1)) )
        self.E_0[0,0] = E_0
        self.F_mat    = F_mat
        self.V_mat    = V_mat
    
    def get_hf_energy(self):
        return self.E_0[0,0]

    def get_fock_mat(self):
        return self.F_mat

    def get_rep_mat(self):
        return self.V_mat

    def __call__(self, cisd_amp_obj):
        return normal_order_hamiltonian_act_on_vec( cisd_amp_obj, self.get_hf_energy(), self.get_fock_mat(), self.get_rep_mat() )




# if __name__ == "__main__":
#     import cisd_amplitude
#     cisd_vec = cisd_amplitude.cisd_amplitude(1,1,4)
#     np.set_printoptions(threshold=np.nan, linewidth=260)
#     #print_amp(cisd_vec)
    
#     def get_fake_F_mat(n):
#         F_mat = np.zeros((2*n,2*n))
#         for i in range(2*n):
#             F_mat[i] = np.random.random(2*n)
#         for i in range(2*n):
#             for j in range(i, 2*n):
#                 F_mat[j,i] = F_mat[i,j]

#         return F_mat

#     def get_fake_V_mat(n):
#         V_mat = np.matrix(np.zeros((4*n*n,4*n*n)))
#         for i in range(4*n*n):
#             V_mat[i] = np.random.random(4*n*n)
#         for i in range(4*n*n):
#             for j in range(i, 4*n*n):
#                 V_mat[j,i] = V_mat[i,j]
#         return V_mat
#     # NEED TO TAKE INPUTS OF E_0, F_{pq} and V_{pqrs}
#     # NEED TO BE FIXED HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     E_0 = -0.77
#     n = cisd_vec.get_num_spatial_orb()
    
#     #F_mat = np.arange(1, 4*n*n+1).reshape(2*n,2*n)
#     # print("Fake Fock Matrix")
#     F_mat = get_fake_F_mat(n)
#     # print(F_mat)

#     # print("Fake V Matrix")
#     V_mat = get_fake_V_mat(n)
#     # print(V_mat)

#     ##################################################################
#     ##################################################################


#     # vec1 = normal_order_hamiltonian_main( cisd_vec, E_0, F_mat, V_mat )
#     # vec2 = normal_order_hamiltonian_main( vec1, E_0, F_mat, V_mat )
#     H = normal_order_hamiltonian( E_0, F_mat, V_mat )
#     # print(type(H.get_hf_energy()))
#     vec_old = cisd_vec
#     i = 0
#     while i < 10:
#         vec_new = H( vec_old )
#         vec_old = vec_new
#         i += 1

#     # print_amp(vec_new)







