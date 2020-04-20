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
import copy



def get_occ_virt_orbs(num_alpha_elec,num_beta_elec,num_spin_orb):
	num_spatial_orb = num_spin_orb // 2
	num_occ_orb  = num_alpha_elec + num_beta_elec
	num_virt_orb = num_spin_orb - num_occ_orb
	spin_orbs    = [ i for i in range(num_spin_orb) ]
	alpha_orbs   = spin_orbs[:num_alpha_elec]
	beta_orbs    = spin_orbs[num_spatial_orb:num_beta_elec + num_spatial_orb]
	alpha_virt_orbs = spin_orbs[num_alpha_elec:num_spatial_orb]
	beta_virt_orbs  = spin_orbs[num_spatial_orb+num_beta_elec:]
	return num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs



def generate_CI_ref_config( num_alpha_elec, num_beta_elec, num_spin_orb ):
	config = [ False for i in range(num_spin_orb) ]
	num_spatial_orb = num_spin_orb // 2
	for i in range(num_alpha_elec):
		config[i] = True
	for i in range(num_beta_elec):
		config[i+num_spatial_orb] = True
	return config



def generate_CI_single_config( num_alpha_elec, num_beta_elec, num_spin_orb ):
	ref_config = generate_CI_ref_config(num_alpha_elec, num_beta_elec, num_spin_orb)
	CIS_config = []
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
			get_occ_virt_orbs(num_alpha_elec,num_beta_elec,num_spin_orb)
	for i in alpha_orbs + beta_orbs:
		for a in alpha_virt_orbs + beta_virt_orbs:
			new_config = copy.deepcopy(ref_config)
			new_config[i] = False
			new_config[a] = True
			CIS_config += [new_config]

	return CIS_config


def generate_CI_double_config( num_alpha_elec, num_beta_elec, num_spin_orb ):
	ref_config = generate_CI_ref_config(num_alpha_elec, num_beta_elec, num_spin_orb)
	CID_config = []
	num_occ_orb, num_virt_orb, alpha_orbs, beta_orbs, alpha_virt_orbs, beta_virt_orbs =\
			get_occ_virt_orbs(num_alpha_elec,num_beta_elec,num_spin_orb)
	occ_orbs  = alpha_orbs + beta_orbs
	virt_orbs = alpha_virt_orbs + beta_virt_orbs
	for i in range(num_occ_orb):
		for j in range(i+1, num_occ_orb):
			for a in range(num_virt_orb):
				for b in range(a+1, num_virt_orb):
					# print(occ_orbs[i], occ_orbs[j],virt_orbs[a], virt_orbs[b])
					new_config = copy.deepcopy(ref_config)
					new_config[ occ_orbs[i] ] = False
					new_config[ occ_orbs[j] ] = False
					new_config[ virt_orbs[a] ] = True
					new_config[ virt_orbs[b] ] = True
					CID_config += [new_config]
	return CID_config


def generate_CISD_configs( num_alpha_elec, num_beta_elec, num_spin_orb, excitation='CISD' ):
	#
	# excitation can be specified as any combinations of the following three,
	# 'CI' ==> Reference
	# 'S'  ==> Singles
	# 'D'  ==> Doubles
	#
	whole_config = []
	if 'CI' in excitation:
		ref_config    = generate_CI_ref_config(num_alpha_elec, num_beta_elec, num_spin_orb)
		whole_config += [ref_config]
	if 'S' in excitation:
		CIS_config    = generate_CI_single_config( num_alpha_elec, num_beta_elec, num_spin_orb )
		whole_config += CIS_config
	if 'D' in excitation:
		CID_config    = generate_CI_double_config( num_alpha_elec, num_beta_elec, num_spin_orb )
		whole_config += CID_config
	return whole_config





if __name__ == "__main__":
	cisd_config = generate_CISD_configs(1,1,8, 'CISD')
	for line in cisd_config:
		print(line)
	print("LENGTH =",len(cisd_config))




