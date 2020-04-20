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
import os
import sys
import numpy
import qode
from qode.util import parallel, output, textlog
import hamiltonian
import ccsd

# Utilities for doing S^(-1/2) orthonormalization
def pos_semi_pow(p):
	def f(x):
		if x<-1e-15:  raise Exception("something horribly wrong?: {}".format(r))
		return abs(x)**p
	return f
def sym_possemi_mpow(M,p):
	Evals,Evecs = numpy.linalg.eigh(M)
	Evals = numpy.vectorize(pos_semi_pow(p))(Evals)
	D = numpy.diag(Evals)
	return Evecs.dot(D.dot(Evecs.T))

# Quick and dirty utility to print to the screen now but append to output later
txt_out = []
def echo_out(string):
	global txt_out
	txt_out += [string]
	print(string)



# Usage:  python3 H.py <basis> <num cation states (1e)>  <num neutral states (2e)>  <num anion states (3e)> <distance> {thresh=<thresh> | nstates=<nstates> | load=<filename>} [hermitian]
basis       = sys.argv[1]
eN          = sys.argv[2:5]
distance    = sys.argv[5]
compression = sys.argv[6].split("=")
hermitian   = (len(sys.argv)==8)

e1, e2, e3      = [int(e) for e in eN]
path            = "/".join(("dimer_H", basis, "-".join(eN)))
out, resources  = output(log=textlog(echo=True)), parallel.resources(1)



echo_out("\n=== info from compression to subspace ===")

# Load/build the non-Hermitian Hamiltonian info for a 2-fragment system and check dimensions
nrg  = numpy.load(path+"/raw/"+distance+"/all_Be_FCI_energies.npy")      	# blocked as |2e>, then |1e>, then |3e> (2e must go first, since the "0" row/column/element is interpreted as the reference in the CC code)
Hmat = numpy.load(path+"/raw/"+distance+"/LMO_fci_H_matrix_raw.npy")		# In non-ON basis. |2e x 2e>, then |2e x 1e>, then |2e x 3e>, then |1e x 2e>, then |1e x 1e>, then |1e x 3e>, then |3e x 2e>, then |3e x 1e>, then |3e x 3e>,
Smat = numpy.load(path+"/raw/"+distance+"/LMO_fci_S_matrix.npy")		# Overlaps.  Same blocking.

for i in range(81):
	row = ""
	for j in range(81):
		if Hmat[i,j]==0:  row += "  "
		else:             row += "X "
	row += "|"
	#echo_out(row)

if hermitian:
	rt_Sinv = sym_possemi_mpow(Smat,-1/2.)
	H = rt_Sinv.dot(Hmat.dot(rt_Sinv))
else:
	Sinv = numpy.linalg.inv(Smat)
	H = Sinv.dot(Hmat)
frag_dim  = len(nrg)
super_dim = frag_dim**2
if frag_dim!=(e1 + e2 + e3):  raise Exception("wrong size energy list loaded for these state numbers")
if H.shape[0]!=super_dim:     raise Exception("wrong size matrix loaded for these state numbers")

# Obtain a small number of states for each monomer
if   compression[0]=="load":

	# Load stored states from disk
	u = numpy.load(compression[1])
	states_per_frag = len(u)

	# Set the path for output (tag references file from with states were located)
	tag = "load=" + compression[1].replace("/",":")
	if hermitian:  tag += "-hermitian"
	out_path = path + "/finished/" + distance + "/" + tag
	if not os.path.exists(out_path):  os.makedirs(out_path)

else:

	# Get lowest dimer energy eigenstate
	Evals,Evecs = numpy.linalg.eig(H)
	Evals,Evecs = zip(*sorted(zip( Evals, Evecs.T.tolist() ), key=lambda p: p[0]))
	Evec, iEvec = numpy.real(Evecs[0]), numpy.imag(Evecs[0])
	if numpy.linalg.norm(iEvec)>1e-15:  raise Exception("something horribly wrong?: norm of eigenvector imaginary part = {}".format(numpy.linalg.norm(iEvec)))
	Evec /= numpy.linalg.norm(Evec)
	echo_out("E_gs  = {}  [from H as loaded]".format(Evals[0]))


	configs = ["|2a>","|2b>","|2c>","|1a>","|1b>","|1c>","|3a>","|3b>","|3c>"]
	Aidx = 0
	Bidx = 0
	for i,c in enumerate(Evec):
		#if   i%9==0:  echo_out("\n==========\n")
		#elif i%3==0:  echo_out(  "----------"  )
		#echo_out("{: 10.2e}   {}{}".format(c,configs[Aidx%9],configs[Bidx%9]))
		Bidx += 1
		if Bidx%9==0:  Aidx += 1

	# Build reduced density matrices for fragments A and B for the dimer ground state
	# [These should be blocked by symmetries, starting with particle number  !!!!!!!!!!! ... ]
	FragA = qode.math.linear_inner_product_space(qode.math.numpy_space.real_traits(frag_dim))
	FragB = qode.math.linear_inner_product_space(qode.math.numpy_space.real_traits(frag_dim))
	Avecs = [ FragB.member(numpy.array(Evec[i:i+frag_dim:1] )) for i in range(0,super_dim,frag_dim) ]
	Bvecs = [ FragA.member(numpy.array(Evec[i:super_dim:frag_dim])) for i in range(0,frag_dim,1)    ]
	rhoA = numpy.array([[(a1|a2) for a2 in Avecs] for a1 in Avecs])
	rhoB = numpy.array([[(b1|b2) for b2 in Bvecs] for b1 in Bvecs])
	rho = (rhoA + rhoB) / 2		# [... !!!!!!!!!!! but this is good enough for linear systems.]

	part_num = [2,2,2, 1,1,1, 3,3,3]

	#echo_out("rho\n")
	for i in range(frag_dim):
		row = ""
		for j in range(frag_dim):
			if part_num[i]!=part_num[j]:  rho[i,j] = 0
			row += "{: 8.0e}".format(rho[i,j])
		#echo_out(row+"\n\n")

	# Diagonalize and sort density matrix eigenvalues/vectors
	rvals, rvecs = numpy.linalg.eigh(rho)
	rvals, rvecs = zip(*reversed(sorted(zip( rvals, rvecs.T.tolist() ), key=lambda p: p[0])))

	#echo_out("natural particle states\n")
	for i in range(frag_dim):
		row = ""
		for j in range(frag_dim):
			row += "{: 8.0e}".format(rvecs[j][i])
		#echo_out(row+"\n\n")

	#utmp = numpy.array(rvecs)
	#rho2 = utmp.dot(rho.dot(utmp.T))
	#echo_out("rebuilt rho\n")
	#for i in range(frag_dim):
	#	row = ""
	#	for j in range(frag_dim):
	#		row += "{: 8.0e}".format(rho2[i,j])
	#	echo_out(row+"\n\n")


	# decide which states to keep
	states_per_frag = 0
	if   compression[0]=="nstates":
		states_per_frag = int(compression[1])
		tag = "n_"+compression[1]
	elif compression[0]=="thresh":
		population_thresh = float(compression[1])
		for r in rvals:
			if r>population_thresh:  states_per_frag += 1
			if r<-1e-15:  raise Exception("something horribly wrong?: {}".format(r))
		tag = compression[1]+"_"+str(states_per_frag)
	else:
		raise Exception("unknown compression algorithm: {}".format(compression[0]))
	echo_out("number of retained fragment states:  {} (out of {})".format(states_per_frag,frag_dim))
	u = list(rvecs[:states_per_frag])

	# Set path for output and dump transformation
	if hermitian:  tag += "-hermitian"
	out_path = path + "/finished/" + distance + "/" + tag
	if not os.path.exists(out_path):  os.makedirs(out_path)
	numpy.save(out_path+"/u.npy", u)


# Make sure the transformations are blocked by particle number with strict zeros (for efficiency, since zeros detected by BCH code)
utmp = numpy.array(u)
ID = utmp.dot(utmp.T)
for i in range(frag_dim):  ID[i,i] -= 1
#print("check!", numpy.linalg.norm(ID))


P1 = numpy.zeros((frag_dim,frag_dim))
P2 = numpy.zeros((frag_dim,frag_dim))
P3 = numpy.zeros((frag_dim,frag_dim))
for i in range(0, e2):            P2[i,i] = 1
for i in range(e2, e2+e1):        P1[i,i] = 1
for i in range(e2+e1, e2+e1+e3):  P3[i,i] = 1
for i in range(len(u)):
	projections = [P.dot(numpy.array(u[i])) for P in (P1,P2,P3)]
	norms = [v.dot(v) for v in projections]
	u[i] = projections[numpy.argmax(norms)]

#u = []
#for i in range(frag_dim):
#	row = [0. for _ in range(frag_dim)]
#	row[i] = 1.
#	u += [row]

utmp = numpy.array(u)
ID = utmp.dot(utmp.T)
for i in range(frag_dim):  ID[i,i] -= 1
#print("check!", numpy.linalg.norm(ID))


# Build the Hamiltonian transformed into the basis of natural monomer states
U = numpy.array([ [a*b for a in uA for b in uB] for uA in u for uB in u ])
ID = U.dot(U.T)
for i in range(super_dim):  ID[i,i] -= 1
#print("check!", numpy.linalg.norm(ID))
h = U.dot(H.dot(U.T))

# Extract the lowest eigenvalue from the compressed Hamiltonian for comparison
Evals,Evecs = numpy.linalg.eig(h)
Evals,Evecs = zip(*sorted(zip( Evals, Evecs.T.tolist() ), key=lambda p: p[0]))
echo_out("E_gs  = {}  [compressed H]".format(Evals[0]))

# Dump compressed Hamiltonian separated into monomer and dimer parts
u,nrg = numpy.array(u),numpy.diag(nrg)
f = u.dot(nrg.dot(u.T))
h0 = numpy.zeros((states_per_frag*states_per_frag, states_per_frag*states_per_frag))
for i in range(states_per_frag):
	for j in range(states_per_frag):
		Mij = i*states_per_frag + j
		for k in range(states_per_frag):
			for l in range(states_per_frag):
				Nkl = k*states_per_frag + l
				if i==k:  h0[Mij][Nkl] += f[j][l]
				if j==l:  h0[Mij][Nkl] += f[i][k]
numpy.save(out_path+"/f.npy", f)
numpy.save(out_path+"/V.npy", h-h0)



# Check our work
out.log("EXCITONIC CCSD CALCULATION")
ccsd.main(hamiltonian.load_subH(path+"/finished", [[".",distance],[".","."]], tag), out, resources)

# Could call fci.main() here, but until everything is rock solid, I like independence of this simple dimer-only H-build
f_mat = numpy.load(out_path+"/f.npy")
V_mat = numpy.load(out_path+"/V.npy")
out.log("Building FCI matrix ...")
for i in range(states_per_frag):
	for j in range(states_per_frag):
		Mij = i*states_per_frag + j
		for k in range(states_per_frag):
			for l in range(states_per_frag):
				Nkl = k*states_per_frag + l
				if i==k:  V_mat[Mij][Nkl] += f_mat[j][l]
				if j==l:  V_mat[Mij][Nkl] += f_mat[i][k]

out.log("Diagonalizing ...")
vals, vecs = numpy.linalg.eig(V_mat)
vals, vecs = zip(*sorted(zip( vals, vecs.T.tolist() ), key=lambda p: p[0]))
out.log("FCI Energy  =", vals[0])

for line in txt_out:  out.log(line)
out.log.write(out_path+"/log.out")
