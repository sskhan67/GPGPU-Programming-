import sys
import pickle
import numpy
from build_density_tensors import build_density_tensors
import the_states

h     = numpy.load("data/semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_h.npy".format(sys.argv[2]))
V_raw = numpy.load("data/semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_V.npy".format(sys.argv[2]))
n_orb = h.shape[0]
V = numpy.zeros((n_orb, n_orb, n_orb, n_orb))
for p in range(n_orb):
    for q in range(n_orb):
        for r in range(n_orb):
            for s in range(n_orb):
                V[p][q][r][s] =  (1/4.) * (V_raw[p][q][r][s] - V_raw[p][q][s][r])



class states(object): pass

neutral = states()
neutral.coeffs  = the_states.neutral
neutral.configs = numpy.array(numpy.load('data/atomic_states/configs_2e.npy'), dtype=numpy.int64)
cation = states()
cation.coeffs  = the_states.cationic
cation.configs = numpy.array(numpy.load('data/atomic_states/configs_1e.npy'), dtype=numpy.int64)
anion = states()
anion.coeffs  = the_states.anionic
anion.configs = numpy.array(numpy.load('data/atomic_states/configs_3e.npy'), dtype=numpy.int64)

z_lists = { 0:neutral , +1:cation, -1:anion }



rho = build_density_tensors(z_lists, n_orbs=9, n_core=1)	# n_orbs and n_core are numbers of spatial orbitals



def monomer(subsys, I, J, rho, h, V):
    n_orb_tot = h.shape[0]	# only for . . . 
    n_orb = n_orb_tot // 2	#  . . . homo-dimeric system
    m, = subsys
    i_m, = I
    j_m, = J
    i_m_chg, i_m_idx = i_m
    j_m_chg, j_m_idx = j_m
    H_IJ = 0
    if i_m_chg != j_m_chg:
        pass
    else:
        Rca   = rho["ca"  ][(i_m_chg, j_m_chg)][i_m_idx][j_m_idx]
        Rccaa = rho["ccaa"][(i_m_chg, j_m_chg)][i_m_idx][j_m_idx]
        H_IJ = 0
        for p in range(n_orb):
            for q in range(n_orb):
                H_IJ += h[m*n_orb+p][m*n_orb+q] * Rca[p][q]
                for r in range(n_orb):
                    for s in range(n_orb):
                        H_IJ += V[m*n_orb+p][m*n_orb+q][m*n_orb+r][m*n_orb+s] * Rccaa[p][q][s][r]
    return H_IJ

def dimer(subsys, I, J, rho, h, V, P2):		# P2 is the parity of the number of electrons in system 2, when it is neutral, 0=even, 1=odd
    n_orb_tot = h.shape[0]	# only for . . . 
    n_orb = n_orb_tot // 2	#  . . . homo-dimeric system
    m1, m2 = subsys
    i_m1, i_m2 = I
    j_m1, j_m2 = J
    i_m1_chg, i_m1_idx = i_m1
    i_m2_chg, i_m2_idx = i_m2
    j_m1_chg, j_m1_idx = j_m1
    j_m2_chg, j_m2_idx = j_m2
    H_IJ = 0
    if i_m1_chg+i_m2_chg != j_m1_chg+j_m2_chg:
        pass
    elif i_m1_chg == j_m1_chg-2:
        Rcc1 = rho["cc"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]	# only for . . .
        Raa2 = rho["aa"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]	#  . . . homo-dimeric system
        for p1 in range(n_orb):
            for q1 in range(n_orb):
                for r2 in range(n_orb):
                    for s2 in range(n_orb):
                        H_IJ += V[m1*n_orb+p1][m1*n_orb+q1][m2*n_orb+r2][m2*n_orb+s2] * Rcc1[p1][q1] * Raa2[s2][r2]
    elif i_m1_chg == j_m1_chg-1:
        Rc1   = rho["c"  ][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]	#
        Rcca1 = rho["cca"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]	# only for . . .
        Ra2   = rho["a"  ][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]	#  . . . homo-dimeric system
        Rcaa2 = rho["caa"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]	#
        sign = -(-1)**((P2+j_m2_chg)%2)
        for p1 in range(n_orb):
            for q2 in range(n_orb):
                H_IJ += sign * h[m1*n_orb+p1][m2*n_orb+q2] * Rc1[p1] * Ra2[q2]
        for p1 in range(n_orb):
            for q1 in range(n_orb):
                for r1 in range(n_orb):
                    for s2 in range(n_orb):
                        H_IJ += sign * 2 * V[m1*n_orb+p1][m1*n_orb+q1][m1*n_orb+r1][m2*n_orb+s2] * Rcca1[q1][p1][r1] * Ra2[s2]
        for p1 in range(n_orb):
            for q2 in range(n_orb):
                for r2 in range(n_orb):
                    for s2 in range(n_orb):
                        H_IJ += sign * 2 * V[m1*n_orb+p1][m2*n_orb+q2][m2*n_orb+r2][m2*n_orb+s2] * Rc1[p1] * Rcaa2[q2][s2][r2]
    elif i_m1_chg == j_m1_chg:
        Rca1 = rho["ca"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]	# only for . . .
        Rca2 = rho["ca"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]	#  . . . homo-dimeric system
        for p1 in range(n_orb):
            for q2 in range(n_orb):
                for r1 in range(n_orb):
                    for s2 in range(n_orb):
                        H_IJ += 4 * V[m1*n_orb+p1][m2*n_orb+q2][m1*n_orb+r1][m2*n_orb+s2] * Rca1[p1][r1] * Rca2[q2][s2]
    elif i_m1_chg == j_m1_chg+1:
        Ra1   = rho["a"  ][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx] 	#
        Rcaa1 = rho["caa"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]	# only for . . .
        Rc2   = rho["c"  ][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]	#  . . . homo-dimeric system
        Rcca2 = rho["cca"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]	#
        sign = +(-1)**((P2+j_m2_chg)%2)
        for p2 in range(n_orb):
            for q1 in range(n_orb):
                H_IJ += sign * h[m2*n_orb+p2][m1*n_orb+q1] * Rc2[p2] * Ra1[q1]
        for p2 in range(n_orb):
            for q2 in range(n_orb):
                for r2 in range(n_orb):
                    for s1 in range(n_orb):
                        H_IJ += sign * 2 * V[m2*n_orb+p2][m2*n_orb+q2][m2*n_orb+r2][m1*n_orb+s1] * Rcca2[q2][p2][r2] * Ra1[s1]
        for p2 in range(n_orb):
            for q1 in range(n_orb):
                for r1 in range(n_orb):
                    for s1 in range(n_orb):
                        H_IJ += sign * 2 * V[m2*n_orb+p2][m1*n_orb+q1][m1*n_orb+r1][m1*n_orb+s1] * Rc2[p2] * Rcaa1[q1][s1][r1]
    elif i_m1_chg == j_m1_chg+2:
        Raa1 = rho["aa"][(i_m1_chg, j_m1_chg)][i_m1_idx][j_m1_idx]	# only for . . .
        Rcc2 = rho["cc"][(i_m2_chg, j_m2_chg)][i_m2_idx][j_m2_idx]	#  . . . homo-dimeric system
        for p2 in range(n_orb):
            for q2 in range(n_orb):
                for r1 in range(n_orb):
                    for s1 in range(n_orb):
                        H_IJ += V[m2*n_orb+p2][m2*n_orb+q2][m1*n_orb+r1][m1*n_orb+s1] * Rcc2[p2][q2] * Raa1[s1][r1]
    else:
        pass
    return H_IJ


bra_chg, ket_chg = sys.argv[1].split("-")
bra_chg0, bra_chg1 = bra_chg.split(",")
if   bra_chg0=="1":
	bra_chg0 = +1
	bra_n0 = z_lists[+1].coeffs.shape[0]
elif bra_chg0=="2":
	bra_chg0 = 0
	bra_n0 = z_lists[0].coeffs.shape[0]
elif bra_chg0=="3":
	bra_chg0 = -1
	bra_n0 = z_lists[-1].coeffs.shape[0]
if   bra_chg1=="1":
	bra_chg1 = +1
	bra_n1 = z_lists[+1].coeffs.shape[0]
elif bra_chg1=="2":
	bra_chg1 = 0
	bra_n1 = z_lists[0].coeffs.shape[0]
elif bra_chg1=="3":
	bra_chg1 = -1
	bra_n1 = z_lists[-1].coeffs.shape[0]
ket_chg0, ket_chg1 = ket_chg.split(",")
if   ket_chg0=="1":
	ket_chg0 = +1
	ket_n0 = z_lists[+1].coeffs.shape[0]
elif ket_chg0=="2":
	ket_chg0 = 0
	ket_n0 = z_lists[0].coeffs.shape[0]
elif ket_chg0=="3":
	ket_chg0 = -1
	ket_n0 = z_lists[-1].coeffs.shape[0]
if   ket_chg1=="1":
	ket_chg1 = +1
	ket_n1 = z_lists[+1].coeffs.shape[0]
elif ket_chg1=="2":
	ket_chg1 = 0
	ket_n1 = z_lists[0].coeffs.shape[0]
elif ket_chg1=="3":
	ket_chg1 = -1
	ket_n1 = z_lists[-1].coeffs.shape[0]

for bra in range(bra_n0):
	string = ""
	for ket in range(ket_n0):
		braHket = monomer((0,), ((bra_chg0,bra),), ((ket_chg0,ket),), rho, h, V)
		string += "{:15.10f}".format(braHket)
	print(string)
print()
for bra in range(bra_n1):
	string = ""
	for ket in range(ket_n1):
		braHket = monomer((1,), ((bra_chg1,bra),), ((ket_chg1,ket),), rho, h, V)
		string += "{:15.10f}".format(braHket)
	print(string)
print()
for bra0 in range(bra_n0):
	for bra1 in range(bra_n1):
		string = ""
		for ket0 in range(ket_n0):
			for ket1 in range(ket_n1):
				braHket = dimer((0,1), ((bra_chg0,bra0),(bra_chg1,bra1)), ((ket_chg0,ket0),(ket_chg1,ket1)), rho, h, V, 0)
				if bra0==ket0:  braHket += monomer((1,), ((bra_chg1,bra1),), ((ket_chg1,ket1),), rho, h, V)
				if bra1==ket1:  braHket += monomer((0,), ((bra_chg0,bra0),), ((ket_chg0,ket0),), rho, h, V)
				string += "{:15.10f}".format(braHket)
		print(string)





