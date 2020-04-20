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
from operator import xor
import copy
import cisd_amplitude
import numpy as np
import qode.SCF_diag_full.hf as scf_routine
import simple_cisd
import cisd_Lanczos
import normal_order_hamiltonian
from math import sqrt
import matrix_diff
# from Applications.general_ci import hamiltonian



def print_amp( cisd_vec ):
    print("=======================================")
    print("-------------\nRef")
    print(cisd_vec.get_ref_amplitude())
    print("-------------\nSingles")
    print(cisd_vec.get_single_amplitude())
    print("-------------\nDoubles")
    print(cisd_vec.get_double_amplitude())
    print("=======================================")

def amp_dot(cisd_vec1, cisd_vec2):
	dot_prod = 0.0
	# Sum over ref
	dot_prod += cisd_vec1.get_ref_amplitude() * cisd_vec2.get_ref_amplitude()
	# 
	# Sum over Singles
	single_mat1 = cisd_vec1.get_single_amplitude()
	single_mat2 = cisd_vec2.get_single_amplitude()
	x_len, y_len = single_mat1.shape
	for i in range(x_len):
		for j in range(y_len):
			dot_prod += single_mat1[i,j] * single_mat2[i,j]
	#
	# Sum over Doubles
	double_mat1 = cisd_vec1.get_double_amplitude()
	double_mat2 = cisd_vec2.get_double_amplitude()
	x_len, y_len = double_mat1.shape
	for i in range(x_len):
		for j in range(i+1, y_len):
			dot_prod += double_mat1[i,j] * double_mat2[i,j]

	return dot_prod


def get_ref_orb_list( num_alpha_elec , num_beta_elec, num_spatial_orb ):
	# return a list of alpha-occ and beta-occ in absolute position number as seen in a alpha + beta dimension.
	ref_list = [ 0  for i in range(num_alpha_elec + num_beta_elec)]
	for i in range(num_alpha_elec):
		ref_list[i] = i
	for i in range(num_beta_elec):
		ref_list[i+num_alpha_elec] = i + num_spatial_orb
	return ref_list


def get_excited_orb_list( ref_list, excitation_list, num_spatial_orb ):
	# list_of_excitations = a 2D list of occ => virt [ [occ1, virt1], [occ2, virt2], ...  ... ]:
	# occ and virt are absolute orbital numbers.
	result_list = copy.deepcopy(ref_list)
	for excitation in excitation_list:
		if excitation[0] < num_spatial_orb:
			result_list[excitation[0]] = excitation[1]
		else:
			result_list[excitation[0] - num_spatial_orb + 1] = excitation[1]
	return result_list


def check_vec_difference( vec1, vec2 ):
	same_orbs = []
	len_vec = len(vec1)
	index_dummy = [i for i in range(len_vec)]
	same_orb_index1 = []
	same_orb_index2 = []
	
	for i in range(len_vec):
		for j in range(len_vec):
			if vec1[i] == vec2[j]:
				same_orbs += [vec1[i]]
				same_orb_index1 += [i]
				same_orb_index2 += [j]
	# print('same orbs = ',same_orbs)
	# print('same in 1 = ',same_orb_index1)
	# print('same in 2 = ',same_orb_index2)
	diff_orb1 = [ i for i in vec1 if i not in same_orbs ]
	diff_orb2 = [ i for i in vec2 if i not in same_orbs ]
	# print("diff 1 = ", diff_orb1)
	# print("diff 2 = ", diff_orb2)
	diff_orb1_index = [ i for i in index_dummy if i not in same_orb_index1 ]
	diff_orb2_index = [ i for i in index_dummy if i not in same_orb_index2 ]
	# print('diff 1 index = ', diff_orb1_index)
	# print('diff 2 index = ', diff_orb2_index)
	return same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2



def check_num_unaligned_index( diff_orb1_index, diff_orb2_index ):
	# Min Diff = 0 which is easy
	# MAX Diff = 2 No higher than that.
	# diff = 1, need to swap the sign
	# diff = 0 or 2, no need to swap the sign.
	common_list = list( set(diff_orb1_index).intersection(diff_orb2_index) )
	return len(diff_orb1_index) - len(common_list)



def index_to_excitation(an_index):
	i = an_index
	if i == 0:
		return "REF"
	elif 1 <= i <= 12:
		return "SINGLE"
	elif 13 <= i <= 27:
		return "DOUBLE"
	else:
		return "ERROR\nERROR\nERROR\n"




# def sum_fock_diag( orbs, F_mat ):
# 	diag_sum = 0.0
# 	for i in orbs:
# 		diag_sum += F_mat[i,i]
# 	return diag_sum


# def get_cisd_amp_from_mat_index( mat_index, cisd_vec_ref_len, cisd_vec_single_len, cisd_vec_double_len, num_alpha_elec, num_beta_elec, num_spatial_orb ):
# 	# if i == 0:  => REF
# 	# elif i < cisd_vec_ref_len + cisd_vec_single_len: => Singles
# 	# else: => Doubles
# 	new_vec = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ).clean_all_amplitude()

# 	if mat_index == 0:
# 		return new_vec.update_ref_amplitude(1.0)
# 	elif mat_index < cisd_vec_ref_len + cisd_vec_single_len:
# 		return new_vec.update_single_amplitude(  )
# 	else:
# 		return "Double"		


#
# Energy calculated from Q-Chem 4.3
#
# SCF   energy in the final basis set = -2.8551604262
# CCSD total energy          = -2.87178917
#
#

E_0, KE, NAE, NRE, ERE, F_mat, T_mat, N_mat, V_mat, vHF, num_alpha_elec, num_beta_elec = scf_routine.main("he.in", 1e-10)


h_mat = T_mat + N_mat

num_spatial_orb = F_mat.shape[0] // 2
# print("NUM SPATIAL ORBS =", num_spatial_orb)

print("KINETIC")
print(T_mat)

print("Nuclear")
print(N_mat)

print("vHF")
print(vHF)

print("Sum of the three")
print(T_mat + N_mat + vHF)


print(E_0)
print(F_mat)

# FOCK_COMPONENT has the following choices:
# 'T', 'N', 'V'
#
#
FOCK_COMPONENT = 'all'
#
#
# RESET EVERYTHING HERE.
#
# (1) T matrix only
#
if FOCK_COMPONENT == 'T':
	E_0   = KE
	h_mat = T_mat
	F_mat = T_mat
	V_mat = np.matrix( np.zeros( V_mat.shape ) )
#
# (2) N matrix only
#
elif FOCK_COMPONENT == 'N':
	E_0   = NAE
	h_mat = N_mat
	F_mat = N_mat
	V_mat = np.matrix( np.zeros( V_mat.shape ) )
#
# (3) V part only
#
elif FOCK_COMPONENT == 'V':
	E_0   = ERE
	h_mat = np.zeros( T_mat.shape )
	F_mat = vHF
	# V mat non-zero
#
# (4) All Components
#
else:
	h_mat = T_mat + N_mat
	F_mat = T_mat + N_mat + vHF
#
#
np.save('F_mat.dat', F_mat)
np.save('V_mat.dat', V_mat)


cisd_amp_vec = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb )
H_cisd       = cisd_Lanczos.build_matrix( E_0, F_mat, V_mat )

# print_amp(H_cisd(cisd_amp_vec))
# 
# print("dot product =", amp_dot(cisd_amp_vec, H_cisd(cisd_amp_vec) ))

# Calculate the total length of CISD vector

num_alpha_virt_orb = cisd_amp_vec.get_num_alpha_virt_orb()
num_beta_virt_orb  = cisd_amp_vec.get_num_bete_virt_orb()
num_spatial_orb    = cisd_amp_vec.get_num_spatial_orb()
num_spin_orb       = 2 * num_spatial_orb
num_occ_orb  = num_alpha_elec + num_beta_elec
num_virt_orb = num_alpha_virt_orb + num_beta_virt_orb


cisd_vec_ref_len = 1

H_1st = []
H_2nd = []


spin_orbitals = [i for i in range(num_spin_orb)]
occ_orbitals  = spin_orbitals[:num_alpha_elec] + spin_orbitals[num_spatial_orb:num_spatial_orb + num_beta_elec]
virt_orbitals = spin_orbitals[num_alpha_elec:num_spatial_orb] + spin_orbitals[num_spatial_orb + num_beta_elec: num_spin_orb]

ref_vec       = get_ref_orb_list(num_alpha_elec,num_beta_elec,num_spatial_orb)

# print("SPIN ORBS =", spin_orbitals)
# print("OCC  = ",occ_orbitals)
# print("VIRT = ",virt_orbitals)

# a_list = get_ref_orb_list(num_alpha_elec,num_beta_elec,num_spatial_orb)
# print(a_list)
# print( get_excited_orb_list(a_list, [[0,5],[1,2]]))

state_config = generate_config.generate_CISD_configs(num_alpha_elec, num_beta_elec, num_spin_orb, 'CISD')
state_obj    = state.state(state_config)


for i_x in range(cisd_vec_ref_len):
	# Loop y starts
	for i_y in range(cisd_vec_ref_len):
		
		#
		#
		# 2nd Quan Starts
		# vec1 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
		# vec2 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
		# vec1.clean_all_amplitude()
		# vec2.clean_all_amplitude()
		# vec1.update_ref_amplitude(1.0)
		# vec2.update_ref_amplitude(1.0)
		# H_2nd += [ amp_dot(vec1, H_cisd(vec2))  ]
		# print("2nd Fisrt H[0,0] = ",H_2nd[0])
		#
		#
		# 1st Quan Starts
		# <R|H|R>
		#
		# H_1st += [ E_0 ]
		#
		#
		h1st_temp = 0.0
		
		same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2 =\
				check_vec_difference(ref_vec, ref_vec)
		# Sum Over Fock
		for i in same_orbs:
			# print("REF-REF, sum over F =", i, 'onto' ,i)
			h1st_temp += h_mat[ i,i ]
		
		
		# Sum Over V
		for i in same_orbs:
			for j in same_orbs:
				# print("REF-REF, sum over V =", i,j,i,j)
				h1st_temp += 0.5 * ( V_mat[ i*num_spin_orb + j, i* num_spin_orb + j ] -\
									 V_mat[ i*num_spin_orb + j, j* num_spin_orb + i] )
		
		H_1st += [ h1st_temp ]
		# print("1st First H[0,0] = ", H_1st[0])
		#
		#
		# End

	# Loop y starts
	for i_y in range(num_occ_orb):
		for j_y in range(num_virt_orb):
			# 2nd Quan Starts
			# vec1 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
			# vec2 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
			# vec1.clean_all_amplitude()
			# vec2.clean_all_amplitude()
			# vec1.update_ref_amplitude(1.0)
			# vec2.single_amp_mat[i_y,j_y] = 1.0
			# H_2nd += [ amp_dot(vec1, H_cisd(vec2)) ]
			#
			#
			# 1st Quan Starts
			# <R|H|S>
			#
			# Differ By One For Sure Not Need To Check
			#
			excitations = [[occ_orbitals[i_y], virt_orbitals[j_y]]]
			single_vec  = get_excited_orb_list( ref_vec, excitations, num_spatial_orb )
			# print("EXCITATIONS = ", excitations)
			# print("SINGLE VEC = ",single_vec)
			same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2 =\
				check_vec_difference(ref_vec, single_vec)
			# num_unaligned_diff_orb = check_num_unaligned_index( diff_orb1_index, diff_orb2_index )
			# print("Differ by ", len(diff_orb1))
			# print("Diff Location = ", diff_orb2_index, diff_orb2_index)
			#
			# 
			# Sum over Fock
			# print("Ref-Single: Fock", diff_orb1[0], diff_orb2[0])
			h1st_temp = h_mat[ diff_orb1[0], diff_orb2[0] ]
			# print("FOCK PART = ", h1st_temp)
			#
			# Sum over V elements
			for i in same_orbs:
				# print("Ref-Single: V SAME ORBS =", diff_orb1[0], i, diff_orb2[0], i)
				h1st_temp += V_mat[diff_orb1[0] * num_spin_orb + i, diff_orb2[0] * num_spin_orb + i] -\
							 V_mat[diff_orb1[0] * num_spin_orb + i, i * num_spin_orb + diff_orb2[0]]
			H_1st += [h1st_temp]
			#
			# End of 1st Quan

	# Loop y starts
	for i_y in range(num_virt_orb):
		for j_y in range(i_y+1, num_virt_orb):
			for k_y in range(num_occ_orb):
				for l_y in range(k_y+1, num_occ_orb):
					# 2nd Quan Starts
					# vec1 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
					# vec2 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
					# vec1.clean_all_amplitude()
					# vec2.clean_all_amplitude()
					# vec1.update_ref_amplitude(1.0)
					# vec2.double_amp_mat[i_y * num_occ_orb + k_y, j_y * num_occ_orb + l_y] = 1.0
					# H_2nd += [ amp_dot(vec1, H_cisd(vec2)) ]
					#
					# 
					# 1st Quan Starts
					# <R|H|D>
					#
					# Differ By Two For Sure.
					excitations = [[occ_orbitals[k_y], virt_orbitals[i_y]], [occ_orbitals[l_y], virt_orbitals[j_y]] ]
					double_vec = get_excited_orb_list( ref_vec, excitations, num_spatial_orb )
					# print("DOUBLE VEC =", double_vec)
					same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2 =\
						check_vec_difference(ref_vec, double_vec)
					# print("SAME ORBS =",  same_orbs)
					# print("REF-DOUBLE sum over V:", diff_orb1, diff_orb2)
					# print("DIFF LOC =", diff_orb1_index, diff_orb2_index)
					# num_diff_loc = check_num_diff_loc( diff_orb1_index, diff_orb2_index )
					H_1st += [  V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[0] * num_spin_orb + diff_orb2[1] ] -\
								V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[1] * num_spin_orb + diff_orb2[0] ] ]
					#
					# End of 1st Quan



for i_x in range(num_occ_orb):
	for j_x in range(num_virt_orb):
		# Loop y starts
		for i_y in range(cisd_vec_ref_len):
			# 2nd Quan Starts
			# vec1 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
			# vec2 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
			# vec1.clean_all_amplitude()
			# vec2.clean_all_amplitude()
			# vec1.single_amp_mat[i_x,j_x] = 1.0
			# vec2.update_ref_amplitude(1.0)
			# H_2nd += [ amp_dot(vec1, H_cisd(vec2)) ]
			#
			#
			# 1st Quan Starts
			# <S|H|R>
			#
			excitations = [[occ_orbitals[i_x], virt_orbitals[j_x]]]
			single_vec  = get_excited_orb_list( ref_vec, excitations, num_spatial_orb )
			# print("EXCITATIONS = ", excitations[0])
			# print("SINGLE VEC = ",single_vec)
			same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2 =\
				check_vec_difference(single_vec, ref_vec)
			# print("Differ by ", len(diff_orb1))
			# print("Diff Location = ", diff_orb1_index, diff_orb2_index)
			#
			#
			# Sum over Fock
			# print("Single-Ref: Fock: ",diff_orb1[0], diff_orb2[0] )
			h1st_temp = h_mat[ diff_orb1[0], diff_orb2[0] ]
			# print("FOCK PART = ", h1st_temp)
			#
			# 
			# Sum over V elements
			for i in same_orbs:
				# print("Single-REF sum over V", diff_orb1[0], i, diff_orb2[0],i )
				h1st_temp += V_mat[ diff_orb1[0] * num_spin_orb + i, diff_orb2[0] * num_spin_orb + i] -\
								V_mat[diff_orb1[0] * num_spin_orb + i, i* num_spin_orb + diff_orb2[0]]
			H_1st += [h1st_temp]
			#
			# End of 1st Quan
			

		# Loop y starts
		for i_y in range(num_occ_orb):
			for j_y in range(num_virt_orb):
				#
				# 2nd Quan Starts
				#
				# vec1 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
				# vec2 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
				# vec1.clean_all_amplitude()
				# vec2.clean_all_amplitude()
				# vec1.single_amp_mat[i_x,j_x] = 1.0
				# vec2.single_amp_mat[i_y,j_y] = 1.0
				# H_2nd += [ amp_dot(vec1, H_cisd(vec2)) ]
				#
				#
				# 1st Quan Starts
				# <S|H|S>
				#
				excitations1 = [ [occ_orbitals[i_x], virt_orbitals[j_x]] ]
				excitations2 = [ [occ_orbitals[i_y], virt_orbitals[j_y]] ]
				single_vec1  = get_excited_orb_list( ref_vec, excitations1, num_spatial_orb )
				single_vec2  = get_excited_orb_list( ref_vec, excitations2, num_spatial_orb )
				# print("SINGLE VECS = ",single_vec1,single_vec2)
				same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2 =\
					check_vec_difference(single_vec1, single_vec2)
				num_unaligned_diff_orb = check_num_unaligned_index( diff_orb1_index, diff_orb2_index )
				# print("UNALIGNED = ", num_unaligned_diff_orb)
				# print("Differ by ", len(diff_orb1))
				# print("Diff Location = ", diff_orb1_index, diff_orb2_index)
				num_diff_orb = len(diff_orb1) 
				h1st_temp = 0.0
				#
				#
				#
				# Same Vector
				if num_diff_orb == 0:
					# print("SAME ORBS = ",same_orbs)
					#
					# Sum over Fock Elements
					for i in same_orbs:
						# print("Single-Single Fock:", i, i)
						h1st_temp += h_mat[ i,i ]
					#
					# Sum over V elements
					for i in same_orbs:
						for j in same_orbs:
							# print("Single-Single Sum over V",i, j, i, j)
							h1st_temp += 0.5 * ( V_mat[ i*num_spin_orb + j, i* num_spin_orb + j ] -\
										V_mat[ i*num_spin_orb + j, j* num_spin_orb + i] )
					H_1st += [ h1st_temp ]
				#
				#
				#
				# One diff orb
				elif num_diff_orb == 1:
					#
					# Sum Over Fock
					# print("Single-Single Fock:", diff_orb1[0], diff_orb2[0] )
					h1st_temp = h_mat[ diff_orb1[0], diff_orb2[0] ]
					#
					# Sum over V
					for i in same_orbs:							
						# print("Single-Single Sum over V",diff_orb1[0], i, diff_orb2[0], i)
						h1st_temp += V_mat[ diff_orb1[0] * num_spin_orb + i, diff_orb2[0] * num_spin_orb + i] -\
										V_mat[diff_orb1[0] * num_spin_orb + i, i* num_spin_orb + diff_orb2[0]]
					#
					#
					# Swap sign if necessary
					if num_unaligned_diff_orb % 2 == 0:
						H_1st += [h1st_temp]
					else:
						H_1st += [-h1st_temp]
				#
				#
				#
				# Two diff orb
				# No need to swap sign at all.
				elif num_diff_orb == 2:
					H_1st += [ V_mat[ diff_orb1[0]*num_spin_orb + diff_orb1[1] , diff_orb2[0] * num_spin_orb + diff_orb2[1] ] -\
										V_mat[ diff_orb1[0]*num_spin_orb + diff_orb1[1] , diff_orb2[1] * num_spin_orb + diff_orb2[0] ] ]
				#
				#
				# Should never go to this branch and this is trivial
				else:
					H_1st += [0.0]
				#
				# End of 1st Quan
				#


		# Loop y starts
		for i_y in range(num_virt_orb):
			for j_y in range(i_y+1, num_virt_orb):
				for k_y in range(num_occ_orb):
					for l_y in range(k_y+1, num_occ_orb):
						# 2nd Quan Starts
						# vec1 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
						# vec2 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
						# vec1.clean_all_amplitude()
						# vec2.clean_all_amplitude()
						# vec1.single_amp_mat[i_x,j_x] = 1.0
						# vec2.double_amp_mat[i_y * num_occ_orb + k_y, j_y * num_occ_orb + l_y] = 1.0
						# H_2nd += [ amp_dot(vec1, H_cisd(vec2)) ]
						#
						#
						# 1st Quan Starts
						# <S|H|D>
						#
						excitations1 = [ [occ_orbitals[i_x], virt_orbitals[j_x]] ]
						excitations2 = [ [occ_orbitals[k_y], virt_orbitals[i_y]], [occ_orbitals[l_y], virt_orbitals[j_y]] ]
						single_vec1  = get_excited_orb_list( ref_vec, excitations1, num_spatial_orb )
						double_vec2  = get_excited_orb_list( ref_vec, excitations2, num_spatial_orb )
						# print("VECS =", single_vec1, double_vec2)
						# Get the diff
						same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2 =\
							check_vec_difference(single_vec1, double_vec2)
						#
						#
						num_diff_orb = len(diff_orb1) 
						num_unaligned_diff_orb = check_num_unaligned_index( diff_orb1_index, diff_orb2_index )
						# print("UNALIGNED = ", num_unaligned_diff_orb)
						# print("NUM OF DIFF ORB =", len(diff_orb1))
						h1st_temp = 0.0
						#
						#
						# Same Vector:
						if num_diff_orb == 0:
							# 
							#
							# Sum Over Fock:
							for i in same_orbs:
								# print("Single-Double Fock:", i, i)
								h1st_temp += h_mat[i,i]
							#
							# Sum Over V elements:
							for i in same_orbs:
								for j in same_orbs:
									# print("Single-Double sum over V", i,j,i,j)
									h1st_temp += 0.5 * ( V_mat[i*num_spin_orb+j, i*num_spin_orb+j] - \
														V_mat[i*num_spin_orb+j, j*num_spin_orb+i]  )
						#
						#
						# Differ by One:
						elif num_diff_orb == 1:
							# print("differ by one")
							#
							# Sum Over Fock:
							h1st_temp = h_mat[ diff_orb1[0], diff_orb2[0] ]
							#
							# Sum Over V elements:
							for i in same_orbs:
								h1st_temp += V_mat[diff_orb1[0]*num_spin_orb + i, diff_orb2[0] * num_spin_orb + i] -\
												V_mat[diff_orb1[0]*num_spin_orb + i, i * num_spin_orb + diff_orb2[0]]
						#
						#
						# Differ by two:
						elif num_diff_orb == 2:
							# print("differ by two")
							h1st_temp = V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[0] * num_spin_orb + diff_orb2[1] ] -\
										 V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[1] * num_spin_orb + diff_orb2[0] ]
						#
						#
						# Differ by more than 2 ==> ZERO.
						else:
							h1st_temp = 0.0
						#
						#
						# Fix the sign according to num_unaligned_diff_orb 
						#
						if num_unaligned_diff_orb % 2 == 0:
							H_1st += [h1st_temp]
						else:
							H_1st += [-h1st_temp]
						#
						#
						# End of 1st Quan


for i_x in range(num_virt_orb):
	for j_x in range(i_x+1, num_virt_orb):
		for k_x in range(num_occ_orb):
			for l_x in range(k_x+1, num_occ_orb):
				# Loop y starts
				for i_y in range(cisd_vec_ref_len):
					# 2nd Quan Starts
					# vec1 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
					# vec2 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
					# vec1.clean_all_amplitude()
					# vec2.clean_all_amplitude()
					# vec1.double_amp_mat[i_x * num_occ_orb + k_x, j_x * num_occ_orb + l_x] = 1.0
					# vec2.update_ref_amplitude(1.0)
					# H_2nd += [ amp_dot(vec1, H_cisd(vec2)) ]
					#
					#
					# 1st Quan Starts
					# <D|H|R>
					# so diff = 2
					#
					#
					excitations = [[occ_orbitals[k_x], virt_orbitals[i_x]], [occ_orbitals[l_x], virt_orbitals[j_x]] ]
					double_vec = get_excited_orb_list( ref_vec, excitations, num_spatial_orb )
					# print("DOUBLE VEC =", double_vec)
					same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2 =\
						check_vec_difference( double_vec, ref_vec )
					# print("SAME ORBS =",  same_orbs)
					# print("DIFF = ", diff_orb1, diff_orb2)
					# print("DIFF LOC =", diff_orb1_index, diff_orb2_index)
					# num_diff_loc = check_num_diff_loc( diff_orb1_index, diff_orb2_index )
					H_1st += [  V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[0] * num_spin_orb + diff_orb2[1] ] -\
								V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[1] * num_spin_orb + diff_orb2[0] ] ]
					#
					#
					# End of 1st Quan


				# Loop y starts
				for i_y in range(num_occ_orb):
					for j_y in range(num_virt_orb):
						# 2nd Quan Starts
						# vec1 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
						# vec2 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
						# vec1.clean_all_amplitude()
						# vec2.clean_all_amplitude()
						# vec1.double_amp_mat[i_x * num_occ_orb + k_x, j_x * num_occ_orb + l_x] = 1.0
						# vec2.single_amp_mat[i_y,j_y] = 1.0
						# H_2nd += [ amp_dot(vec1, H_cisd(vec2)) ]
						#
						#
						#
						# 1st Quan Starts
						# <D|H|S>
						#
						excitations1 = [ [occ_orbitals[k_x], virt_orbitals[i_x]], [occ_orbitals[l_x], virt_orbitals[j_x]] ]
						excitations2 = [ [occ_orbitals[i_y], virt_orbitals[j_y]] ]
						double_vec1  = get_excited_orb_list( ref_vec, excitations1, num_spatial_orb )
						single_vec2  = get_excited_orb_list( ref_vec, excitations2, num_spatial_orb )
						# print("VECS =", single_vec1, double_vec2)
						# Get the diff
						same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2 =\
							check_vec_difference(double_vec1, single_vec2)
						#
						#
						num_diff_orb = len(diff_orb1) 
						num_unaligned_diff_orb = check_num_unaligned_index( diff_orb1_index, diff_orb2_index )
						# print("UNALIGNED = ", num_unaligned_diff_orb)
						# print("NUM OF DIFF ORB =", len(diff_orb1))
						h1st_temp = 0.0
						#
						#
						# Single and Double are NEVER gonna be the same.
						#
						#
						# Differ by One:
						if num_diff_orb == 1:
							# print("differ by one")
							#
							# Sum Over Fock:
							h1st_temp = h_mat[diff_orb1[0], diff_orb2[0]]
							#
							# Sum Over V elements:
							for i in same_orbs:
								h1st_temp += V_mat[diff_orb1[0]*num_spin_orb + i, diff_orb2[0] * num_spin_orb + i] -\
											 V_mat[diff_orb1[0]*num_spin_orb + i, i * num_spin_orb + diff_orb2[0]]
						#
						#
						# Differ by two:
						elif num_diff_orb == 2:
							# print("differ by two")
							h1st_temp = V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[0] * num_spin_orb + diff_orb2[1] ] -\
										V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[1] * num_spin_orb + diff_orb2[0] ]
						#
						#
						# Differ by more than 2 ==> ZERO.
						# Single and double are never identical
						else:
							h1st_temp = 0.0
						#
						#
						# Fix the sign according to num_unaligned_diff_orb 
						#
						if num_unaligned_diff_orb % 2 == 0:
							H_1st += [h1st_temp]
						else:
							H_1st += [-h1st_temp]
						#
						#
						# End of 1st Quan





				# Loop y starts
				for i_y in range(num_virt_orb):
					for j_y in range(i_y+1, num_virt_orb):
						for k_y in range(num_occ_orb):
							for l_y in range(k_y+1, num_occ_orb):
								# 2nd Quan Starts
								# vec1 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
								# vec2 = cisd_Lanczos.build_vector( num_alpha_elec, num_beta_elec, num_spatial_orb ) 
								# vec1.clean_all_amplitude()
								# vec2.clean_all_amplitude()
								# vec1.double_amp_mat[i_x * num_occ_orb + k_x, j_x * num_occ_orb + l_x] = 1.0
								# vec2.double_amp_mat[i_y * num_occ_orb + k_y, j_y * num_occ_orb + l_y] = 1.0
								# H_2nd += [ amp_dot(vec1, H_cisd(vec2)) ]
								#
								#
								# 1st Quan Starts
								# <D|H|D>
								# 
								# 
								#
								excitations1 = [ [occ_orbitals[k_x], virt_orbitals[i_x]], [occ_orbitals[l_x], virt_orbitals[j_x]] ]
								excitations2 = [ [occ_orbitals[k_y], virt_orbitals[i_y]], [occ_orbitals[l_y], virt_orbitals[j_y]] ]
								double_vec1  = get_excited_orb_list( ref_vec, excitations1, num_spatial_orb )
								double_vec2  = get_excited_orb_list( ref_vec, excitations2, num_spatial_orb )
								# print("VECS =", single_vec1, double_vec2)
								# Get the diff
								same_orbs, diff_orb1_index, diff_orb1 ,diff_orb2_index, diff_orb2 =\
									check_vec_difference(double_vec1, double_vec2)
								#
								#
								num_diff_orb = len(diff_orb1) 
								num_unaligned_diff_orb = check_num_unaligned_index( diff_orb1_index, diff_orb2_index )
								# print("UNALIGNED = ", num_unaligned_diff_orb)
								# print("NUM OF DIFF ORB =", len(diff_orb1))
								h1st_temp = 0.0
								#
								#
								#
								# Same Vector:
								if num_diff_orb == 0:
									#
									# Sum Over Fock:
									for i in same_orbs:
										h1st_temp += h_mat[i,i]
									#
									# Sum Over V elements:
									for i in same_orbs:
										for j in same_orbs:
											h1st_temp += 0.5 * ( V_mat[i*num_spin_orb+j, i*num_spin_orb+j] - \
																 V_mat[i*num_spin_orb+j, j*num_spin_orb+i]  )
								#
								#
								# Differ by One:
								elif num_diff_orb == 1:
									#
									# Sum Over Fock:
									h1st_temp = h_mat[diff_orb1[0], diff_orb2[0]]
									#
									# Sum Over V elements:
									for i in same_orbs:
										h1st_temp += V_mat[diff_orb1[0]*num_spin_orb + i, diff_orb2[0] * num_spin_orb + i] -\
													 V_mat[diff_orb1[0]*num_spin_orb + i, i * num_spin_orb + diff_orb2[0]]
								#
								#
								# Differ by Two:
								elif num_diff_orb == 2:
									h1st_temp = V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[0] * num_spin_orb + diff_orb2[1] ] -\
										 		V_mat[ diff_orb1[0] * num_spin_orb + diff_orb1[1], diff_orb2[1] * num_spin_orb + diff_orb2[0] ]
								#
								#
								# Differ more than two => ZERO
								else:
									h1st_temp = 0.0
								#
								#
								# Fix the sign according to num_unaligned_diff_orb 
								#
								if num_unaligned_diff_orb % 2 == 0:
									H_1st += [h1st_temp]
								else:
									H_1st += [-h1st_temp]
								#
								#
								# End of 1st Quan



# print(H_2nd)
# H_2nd = np.matrix( H_2nd )
# print(H_2nd)
# size2 = H_2nd.shape[1]
# size = int(sqrt( size2 ))
# H_2nd = np.reshape(H_2nd, (size,size))
# E,V = np.linalg.eigh(H_2nd)
# print("2nd Quantization Lowest E = ",E[0])
# print("First Numpy Vec = ", V[0])
# print("H made by 1st QM")
# print(H_1st)
size  = int(sqrt(len(H_1st)))
H_1st = np.array(H_1st)
H_1st = np.reshape(H_1st,(size,size))
H_1st = np.matrix(H_1st)
E1st, V1st = np.linalg.eigh(H_1st)
print("1st Quantization Lowest E = ", E1st[0])
print("CCSD total energy         =  -2.87178917")



np.save("first_quantized_mat.npy", H_1st)


# output = ""
# for i in range(size):
# 	for j in range(size):
# 		i_state = index_to_excitation(i)
# 		j_state = index_to_excitation(j)
# 		if abs(H_2nd[i,j] - H_1st[i,j]) > 1e-15:
# 			output += "H[%d,%d]  %s %s First = %e  Second = %e Diff = %e"  %( i,j,i_state, j_state ,H_1st[i,j] ,H_2nd[i,j], H_2nd[i,j] - H_1st[i,j] )
# 			if abs(H_2nd[i,j]) == abs(H_1st[i,j]):
# 				output += "  OPPOSITE SIGN\n"
# 			else:
# 				output += '\n'
# 		# else:
# 		# 	output += '\n'
# f = open(FOCK_COMPONENT + "_diff_mat.txt",'w')
# f.write(output)
# f.close()
# print("FILE WRITTEN TO diff_mat.txt")

# # matrix_diff.matrix_sign_diff(H_1st, H_2nd)
# print( matrix_diff.cmp_mat(H_1st,H_2nd) )

np.set_printoptions(precision=3,linewidth=270,threshold=np.nan)
print("FIRST Quantization")
print(H_1st)
# print("Second Quantization")
# print(H_2nd)









