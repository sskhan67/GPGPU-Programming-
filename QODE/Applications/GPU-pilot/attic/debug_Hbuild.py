import sys
import numpy
from qode.fermion_field.state import *

h1     = numpy.load("data/semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_h.npy".format(sys.argv[2]))		# <^p| h |_q>
v2_raw = numpy.load("data/semiMO_ints/{}_A/biorthogonal_semiMO_631G_Be2_V.npy".format(sys.argv[2]))		# <^p| <^q| V |_r> |_s>
n_orb = h1.shape[0]
v2 = numpy.zeros((n_orb, n_orb, n_orb, n_orb))
for p in range(n_orb):
	for q in range(n_orb):
		for r in range(n_orb):
			for s in range(n_orb):
				v2[p][q][r][s] =  v2_raw[p][q][r][s] - v2_raw[p][q][s][r]			# <^p^q| V |_r_s>



nsys, bra_chg, ket_chg = sys.argv[1].split("-")
if   nsys=="monomer":
	if bra_chg!=ket_chg:  raise
	orb_lists = numpy.load("data/atomic_states/configs_{}e.npy".format(bra_chg))
	occ_lists = []
	for orb_list in orb_lists:
		occ_list = [False]*18
		for orb in orb_list:  occ_list[orb] = True
		occ_lists += [occ_list]
	BRA_OCC_LISTS = occ_lists
	KET_OCC_LISTS = occ_lists
elif nsys=="dimer":
	bra_chgA, bra_chgB = bra_chg.split(",")
	ket_chgA, ket_chgB = ket_chg.split(",")
	if int(bra_chgA)+int(bra_chgB) != int(ket_chgA)+int(ket_chgB):  raise
	orb_lists = numpy.load("data/atomic_states/configs_{}e.npy".format(bra_chgA))
	occ_listsA = []
	for orb_list in orb_lists:
		occ_list = [False]*18
		for orb in orb_list:  occ_list[orb] = True
		occ_listsA += [occ_list]
	orb_lists = numpy.load("data/atomic_states/configs_{}e.npy".format(bra_chgB))
	occ_listsB = []
	for orb_list in orb_lists:
		occ_list = [False]*18
		for orb in orb_list:  occ_list[orb] = True
		occ_listsB += [occ_list]
	BRA_OCC_LISTS = []
	for occ_listA in occ_listsA:
		for occ_listB in occ_listsB:
			BRA_OCC_LISTS += [ occ_listA + occ_listB ]
	orb_lists = numpy.load("data/atomic_states/configs_{}e.npy".format(ket_chgA))
	occ_listsA = []
	for orb_list in orb_lists:
		occ_list = [False]*18
		for orb in orb_list:  occ_list[orb] = True
		occ_listsA += [occ_list]
	orb_lists = numpy.load("data/atomic_states/configs_{}e.npy".format(ket_chgB))
	occ_listsB = []
	for orb_list in orb_lists:
		occ_list = [False]*18
		for orb in orb_list:  occ_list[orb] = True
		occ_listsB += [occ_list]
	KET_OCC_LISTS = []
	for occ_listA in occ_listsA:
		for occ_listB in occ_listsB:
			KET_OCC_LISTS += [ occ_listA + occ_listB ]
else:  raise



def H_element_0(bra_orbs, ket_orbs, bra, ket):
	Hmn = 0
	for p in bra_orbs:
		for q in ket_orbs:
			op = op_string(create(p),annihilate(q))
			pq_ket = op | ket
			Hmn += (bra|pq_ket) * h1[p,q]
	for p in bra_orbs:
		for q in bra_orbs:
			for r in ket_orbs:
				for s in ket_orbs:
					op = op_string(create(p),create(q),annihilate(s),annihilate(r))
					pqsr_ket = op | ket
					Hmn += (bra|pqsr_ket) * v2[p,q,r,s]/4
	return Hmn

def H_element_1(p, q, bra_orbs, ket_orbs, bra, ket):
	op = op_string(create(p),annihilate(q))
	pq_ket = op | ket
	Hmn = (bra|pq_ket) * h1[p,q]
	for r in bra_orbs:
		op = op_string(create(p),create(r),annihilate(r),annihilate(q))
		prrq_ket = op | ket
		Hmn += (bra|prrq_ket) * v2[p,r,q,r]
	return Hmn

def H_element_2(p, q, r, s, bra, ket):
	op = op_string(create(p),create(q),annihilate(s),annihilate(r))
	pqsr_ket = op | ket
	Hmn = (bra|pqsr_ket) * v2[p,q,r,s]
	return Hmn



H = numpy.zeros((len(BRA_OCC_LISTS),len(KET_OCC_LISTS)))

for m,bra_occ in enumerate(BRA_OCC_LISTS):
	for n,ket_occ in enumerate(KET_OCC_LISTS):
		bra_orbs = []
		ket_orbs = []
		bra_XS   = []
		ket_XS   = []
		diffs    = 0
		for i,(b,k) in enumerate(zip(bra_occ,ket_occ)):
			if b: bra_orbs += [i]	
			if k: ket_orbs += [i]
			if b!=k:
				if b: bra_XS += [i]	
				if k: ket_XS += [i]
				diffs += 1
			if diffs>4:  break
		if   diffs==0:
			bra = configuration(bra_occ)
			ket = configuration(ket_occ)
			H[m,n] = H_element_0(bra_orbs, ket_orbs, bra, ket)
		elif diffs==2:
			p, = bra_XS
			q, = ket_XS
			bra = configuration(bra_occ)
			ket = configuration(ket_occ)
			H[m,n] = H_element_1(p, q, bra_orbs, ket_orbs, bra, ket)
		elif diffs==4:
			p,q = bra_XS
			r,s = ket_XS
			bra = configuration(bra_occ)
			ket = configuration(ket_occ)
			H[m,n] = H_element_2(p, q, r, s, bra, ket)

numpy.save("H/{}_{}.npy".format(sys.argv[1],sys.argv[2]),H)
