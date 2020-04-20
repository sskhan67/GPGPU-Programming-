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
# An object to hold CISD amplitudes, the base class is to be inherited by children.
import sys
import numpy as np
import copy

class amplitude:
    pass


class cisd_amplitude(amplitude):
    def __init__(self,
                 num_alpha_elec,
                 num_beta_elec,
                 num_basis_func
                 ):
        # num_basis_func is the number of spatial orbitals
        self.num_alpha_elec = num_alpha_elec
        self.num_beta_elec  = num_beta_elec
        self.num_basis_func = num_basis_func
        self.num_alpha_virt_orb = self.num_basis_func - self.num_alpha_elec
        self.num_beta_virt_orb  = self.num_basis_func - self.num_beta_elec
        
        # CHECK FIRST IF ALL KINDS OF DOUBLE EXCITATIONS ARE POSSIBLE!!
        if self.num_alpha_virt_orb >= 0 and self.num_beta_virt_orb >= 0:
            pass
        else:
            print("Something really weird went wrong")
            sys.exit(1)
                
        # build three amplitude holders, Ref is for one number,
        # Singles is a matrix completely filled with numbers,
        # Doubles is a list of upper submatrices ( Please refer to the png file in current folder )
        
        # Reference Amplitude:
        self.ref_amplitude = np.matrix(np.zeros((1,1)))
        self.ref_amplitude[0,0] = 1.0
    
        # Singles Amplitudes: ( ROWS: m_alpha , m_beta ; COLUMNS: e_alpha , e_beta )
        self.single_amp_mat = np.matrix(np.zeros(( self.num_alpha_elec + self.num_beta_elec , self.num_alpha_virt_orb \
                                        + self.num_beta_virt_orb )) )
        
        # Doubles Amplitudes: ( ROWS: e_alpha, e_beta; COLUMNS: f_alpha, f_beta )
        # Just make a super buge matrix to hold zeros.
        total_ld = ( self.num_alpha_elec + self.num_beta_elec ) * ( self.num_alpha_virt_orb + self.num_beta_virt_orb )
        self.double_amp_mat = np.matrix( np.zeros(( total_ld, total_ld )) )
    
    # Functions to retrieve electronic configurations
    def get_num_alpha_elec(self):
        return self.num_alpha_elec
    
    def get_num_beta_elec(self):
        return self.num_beta_elec
    
    def get_num_spatial_orb(self):
        return self.num_basis_func
    
    def get_num_alpha_virt_orb(self):
        return self.num_alpha_virt_orb
    
    def get_num_bete_virt_orb(self):
        return self.num_beta_virt_orb
    
    # Functions to retrieve CISD amplitudes
    def get_ref_amplitude(self):
        return self.ref_amplitude[0,0]
           
    def get_single_amplitude(self):
        return self.single_amp_mat
    
    def get_double_amplitude(self):
        return self.double_amp_mat

    # Functions to update CISD amplitudes
    def update_ref_amplitude(self, a_number ):
        self.ref_amplitude[0,0] = a_number
        
    def update_single_amplitude(self, a_np_matrix ):
        if self.single_amp_mat.shape == a_np_matrix.shape:
            self.single_amp_mat = np.matrix( copy.deepcopy( a_np_matrix ) )
        else:
            print("Incompatible Matrix Size assigned to Single Excitation Amplitude Matrix.")
            sys.exit(1)
    
    def update_double_amplitude(self, a_np_matrix ):
        if self.double_amp_mat.shape == a_np_matrix.shape:
            self.double_amp_mat = np.matrix( copy.deepcopy( a_np_matrix ) )
        else:
            print("Incompatible Matrix Size assigned to Double Excitation Amplitude Matrix.")
            sys.exit(1)
    def get_vec_dimension(self):
        return 1 + self.get_single_amplitude().shape[0] * self.get_single_amplitude().shape[1] + \
            self.get_double_amplitude().shape[0] *  self.get_double_amplitude().shape[1]


    def clean_all_amplitude(self):
        self.ref_amplitude[0,0] = 0.0
        self.single_amp_mat.fill(0.0)
        self.double_amp_mat.fill(0.0)



if __name__ == "__main__":
    cisd_vec = cisd_amplitude(1,1,4)
    np.set_printoptions(threshold=np.nan)
    print("-------------\nalpha Ref")
    print(cisd_vec.get_ref_amplitude())
    print("-------------\nalpha singles")
    print(cisd_vec.get_single_amplitude())
    print("-------------\nalpha doubles")
    print(cisd_vec.get_double_amplitude())
    print("TYPE =", type(cisd_vec.get_ref_amplitude()))
    
    