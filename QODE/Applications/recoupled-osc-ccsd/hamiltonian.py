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
from qode.coupled_oscillators.CI_util.generate_ho_config import diag_rec_v_matrix



def load(H, states_per_frag, all_eigenvalues, v_mat_list):

	# 1-body terms
	for n,evals in enumerate(all_eigenvalues):
		# Energy of the reference
		H.K[0] += evals[0]
		# Fo Parts:  <0|h|0>
		H.Fo[n] = -evals[0]
		# Fv Parts: <a|h|a>
		for a,ev in enumerate(evals[1:]):  H.Fv[ H.Fv_starters[n] + a*len(evals) ] += ev
		# Ex Parts: All Zero
		# Dx Parts: All Zero

	# 2-body terms
	num_mol = len(states_per_frag)
	Nv = [ num-1 for num in states_per_frag ]
	V_mat = diag_rec_v_matrix(states_per_frag, v_mat_list)
	#
	for m in range(num_mol):
		for n in range(num_mol):
			if m > n:
				H.K[0]  += V_mat( n, m, 0, 0, 0, 0 )
				H.Fo[m] += -V_mat( n, m, 0, 0, 0, 0 )
			if m < n: 
				H.Fo[m] += -V_mat( m, n, 0, 0, 0, 0 )
	for m in range(num_mol):
		for n in range(num_mol): 
			if m > n:
				for a in range(Nv[m]):
					for b in range(Nv[m]):
						H.Fv[ H.Fv_starters[m] + b*Nv[m] + a ] += V_mat( n, m, 0, a+1, 0, b+1 )
			if m < n:
				for a in range(Nv[m]):
					for b in range(Nv[m]):
						H.Fv[ H.Fv_starters[m] + b*Nv[m] + a ] += V_mat( m, n, a+1, 0, b+1, 0 )
	for m in range(num_mol):
		for n in range(num_mol):
			if m > n: 
				for a in range(Nv[m]):
					H.Ex[ H.Ex_starters[m] + a ] += V_mat( n, m, 0, 0, 0, a+1 )
			if m < n: 
				for a in range(Nv[m]):
					H.Ex[ H.Ex_starters[m] + a ] += V_mat( m, n, 0, 0, a+1, 0 )
	for m in range(num_mol):
		for n in range(num_mol):
			if m > n: 
				for a in range(Nv[m]):
					H.Dx[ H.Dx_starters[m] + a ] += -V_mat( n, m, 0, a+1, 0, 0 )
			if m < n: 
				for a in range(Nv[m]):
					H.Dx[ H.Dx_starters[m] + a ] += -V_mat( m, n, a+1, 0, 0, 0 )
	#
	# ExEx[i,a,j,b]  m<n: <ij|V|ab>; m>n: <ji|V|ba>
	for m in range(num_mol):
		for n in range(num_mol): 
			if m > n: 
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						H.ExEx[ H.ExEx_starters[m*num_mol+n] + a*Nv[n] + b ] += V_mat( n, m, 0, 0, b+1, a+1 )
						H.ExEx[ H.ExEx_starters[n*num_mol+m] + b*Nv[m] + a ] += V_mat( n, m, 0, 0, b+1, a+1 )
	# ExFo: ExFo[i,a,j,k]  m<n: <ij|V|ak>; m>n: <ji|V|ka>
	for m in range(num_mol):
		for n in range(num_mol): 
			if m < n:
				for a in range(Nv[m]):
					H.ExFo[ H.ExFo_starters[m*num_mol+n] + a ] += -V_mat( m, n, 0, 0, a+1, 0 )
			elif m > n:
				for a in range(Nv[m]):
					H.ExFo[ H.ExFo_starters[m*num_mol+n] + a ] += -V_mat( n, m, 0, 0, 0, a+1 )
	# ExFv: ExFv[i,a,b,c]  m<n: <ib|V|ac>;  m>n: <bi|V|ca>
	for m in range(num_mol):
		for n in range(num_mol): 
			if m < n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						for c in range(Nv[n]):
							H.ExFv[ H.ExFv_starters[m*num_mol+n] + a*Nv[n]**2 + c*Nv[n] + b ] += V_mat( m, n, 0, b+1, a+1, c+1 )
			if m > n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						for c in range(Nv[n]):
							H.ExFv[ H.ExFv_starters[m*num_mol+n] + a*Nv[n]**2 + c*Nv[n] + b ] += V_mat( n, m, b+1, 0, c+1, a+1 )
	# ExDx ExDx[i,a,b,j]  <ib|V|aj>   <bi|V|ja>
	for m in range(num_mol):
		for n in range(num_mol):
			if m < n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						H.ExDx[ H.ExDx_starters[m*num_mol+n] + a*Nv[n] + b ] += -V_mat( m, n, 0, b+1, a+1, 0 )
			elif m > n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						H.ExDx[ H.ExDx_starters[m*num_mol+n] + a*Nv[n] + b ] += -V_mat( n, m, b+1, 0, 0, a+1 )
	# FoFo: FoFo[i,j,k,l]  <ik|V|jl>   <ki|V|lj>
	for m in range(num_mol):
		for n in range(num_mol):
			if m > n:
				H.FoFo[ H.FoFo_starters[m*num_mol+n] ] += V_mat( n, m, 0, 0, 0, 0 )
				H.FoFo[ H.FoFo_starters[n*num_mol+m] ] += V_mat( n, m, 0, 0, 0, 0 )
	# FoFv:  FoFv[i,j,a,b]  <ia|V|jb>  <ai|V|bj>
	for m in range(num_mol):
		for n in range(num_mol):
			if m < n:
				for a in range(Nv[n]):
					for b in range(Nv[n]):
						H.FoFv[ H.FoFv_starters[m*num_mol+n] + b*Nv[n] + a ] += -V_mat( m, n, 0, a+1, 0, b+1 )
			elif m > n:
				for a in range(Nv[n]):
					for b in range(Nv[n]):
						H.FoFv[ H.FoFv_starters[m*num_mol+n] + b*Nv[n] + a ] += -V_mat( n, m, a+1, 0, b+1, 0 )
	# FoDx: FoDx[i,j,a, k]   <ia|V|jk>  <ai|V|kj>
	for m in range(num_mol):
		for n in range(num_mol):
			if m < n:
				for a in range(Nv[n]):
					H.FoDx[ H.FoDx_starters[m*num_mol+n] + a ] += V_mat( m, n, 0, a+1, 0, 0 )
			elif m > n:
				for a in range(Nv[n]):
					H.FoDx[ H.FoDx_starters[m*num_mol+n] + a ] += V_mat( n, m, a+1, 0, 0, 0 )
	# FvFv: FvFv[a,b,c,d]   <ac|V|bd>   <ca|V|db>
	for m in range(num_mol):
		for n in range(num_mol): 				
			if m > n:
				for a in range(Nv[m]):
					for b in range(Nv[m]):
						for c in range(Nv[n]):
							for d in range(Nv[n]):
								H.FvFv[ H.FvFv_starters[m*num_mol+n] +  b*Nv[m]*Nv[n]**2 + a*Nv[n]**2 + d*Nv[n] + c ] += V_mat( n, m, c+1, a+1, d+1, b+1 )
								H.FvFv[ H.FvFv_starters[n*num_mol+m] +  d*Nv[n]*Nv[m]**2 + c*Nv[m]**2 + b*Nv[m] + a ] += V_mat( n, m, c+1, a+1, d+1, b+1 )
	# FvDx: FvDx[a,b,c,i]  <ac|V|bi>  <ca|V|ib>
	for m in range(num_mol):
		for n in range(num_mol): 		
			if m < n:		
				for a in range(Nv[m]):
					for b in range(Nv[m]):
						for c in range(Nv[n]):
							H.FvDx[ H.FvDx_starters[m*num_mol+n] + b*Nv[m]*Nv[n] + a*Nv[n] + c ] += -V_mat( m, n, a+1, c+1, b+1, 0 )
			elif m > n:
				for a in range(Nv[m]):
					for b in range(Nv[m]):
						for c in range(Nv[n]):
							H.FvDx[ H.FvDx_starters[m*num_mol+n] + b*Nv[m]*Nv[n] + a*Nv[n] + c ] += -V_mat( n, m, c+1, a+1, 0, b+1 )
	# DxDx: DxDx[a,i,b,j]  <ab|V|ij>  <ba|V|ji>
	for m in range(num_mol):
		for n in range(num_mol):
			if m > n:
				for a in range(Nv[m]):
					for b in range(Nv[n]):
						H.DxDx[ H.DxDx_starters[m*num_mol+n] + a*Nv[n] + b ] += V_mat( n, m, b+1, a+1, 0, 0 )
						H.DxDx[ H.DxDx_starters[n*num_mol+m] + b*Nv[m] + a ] += V_mat( n, m, b+1, a+1, 0, 0 )
