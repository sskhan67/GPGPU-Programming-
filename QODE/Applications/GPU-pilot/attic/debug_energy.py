import sys
import numpy
import the_states

coeffs = { "1":the_states.cationic, "2":the_states.neutral, "3":the_states.anionic }

nsys, bra_chg, ket_chg = sys.argv[1].split("-")

if   nsys=="monomer":
	if bra_chg!=ket_chg:  raise
	BRA_COEFFS = coeffs[bra_chg]
	KET_COEFFS = coeffs[ket_chg]
elif nsys=="dimer":
	bra_chgA, bra_chgB = bra_chg.split(",")
	ket_chgA, ket_chgB = ket_chg.split(",")
	if int(bra_chgA)+int(bra_chgB) != int(ket_chgA)+int(ket_chgB):  raise
	BRA_COEFFS = []
	for coeffsA in coeffs[bra_chgA]:
		for coeffsB in coeffs[bra_chgB]:
			coeffsAB = []
			for coeffA in coeffsA:
				for coeffB in coeffsB:
					coeffsAB += [ coeffA * coeffB ]
			BRA_COEFFS += [coeffsAB]
	KET_COEFFS = []
	for coeffsA in coeffs[ket_chgA]:
		for coeffsB in coeffs[ket_chgB]:
			coeffsAB = []
			for coeffA in coeffsA:
				for coeffB in coeffsB:
					coeffsAB += [ coeffA * coeffB ]
			KET_COEFFS += [coeffsAB]
else:  raise
BRA_COEFFS = numpy.array(BRA_COEFFS)
KET_COEFFS = numpy.array(KET_COEFFS)

H = numpy.load("H/{}_{}.npy".format(sys.argv[1],sys.argv[2]))

for bra in BRA_COEFFS:
	string = ""
	for ket in KET_COEFFS:
		Hket = H.dot(ket)
		braHket = bra.dot(Hket)
		string += "{:15.10f}".format(braHket)
	print(string)
