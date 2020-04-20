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
from data.extract import X
from data import atom

def parse(n, fields):
	subfields = float(fields[0].split("_")[3])
	energy    = float(fields[-1])
	return subfields, energy-2*atom.energy

data = X(
parse,
"""\
be2_fci_6-31g_3.000_A.out:* ROOT 1 CI total energy = -29.2222388231744
be2_fci_6-31g_3.100_A.out:* ROOT 1 CI total energy = -29.2228041600158
be2_fci_6-31g_3.200_A.out:* ROOT 1 CI total energy = -29.2233213238590
be2_fci_6-31g_3.300_A.out:* ROOT 1 CI total energy = -29.2237853390280
be2_fci_6-31g_3.400_A.out:* ROOT 1 CI total energy = -29.2241931867335
be2_fci_6-31g_3.500_A.out:* ROOT 1 CI total energy = -29.2245443827129
be2_fci_6-31g_3.600_A.out:* ROOT 1 CI total energy = -29.2248407498656
be2_fci_6-31g_3.700_A.out:* ROOT 1 CI total energy = -29.2250858632411
be2_fci_6-31g_3.800_A.out:* ROOT 1 CI total energy = -29.2252844337291
be2_fci_6-31g_3.900_A.out:* ROOT 1 CI total energy = -29.2254417597904
be2_fci_6-31g_4.000_A.out:* ROOT 1 CI total energy = -29.2255632959796
be2_fci_6-31g_4.100_A.out:* ROOT 1 CI total energy = -29.2256543440307
be2_fci_6-31g_4.200_A.out:* ROOT 1 CI total energy = -29.2257198521794
be2_fci_6-31g_4.300_A.out:* ROOT 1 CI total energy = -29.2257643011786
be2_fci_6-31g_4.400_A.out:* ROOT 1 CI total energy = -29.2257916549210
be2_fci_6-31g_4.500_A.out:* ROOT 1 CI total energy = -29.2258053560890
be2_fci_6-31g_4.600_A.out:* ROOT 1 CI total energy = -29.2258083508322
be2_fci_6-31g_4.700_A.out:* ROOT 1 CI total energy = -29.2258031301391
be2_fci_6-31g_4.800_A.out:* ROOT 1 CI total energy = -29.2257917788818
be2_fci_6-31g_4.900_A.out:* ROOT 1 CI total energy = -29.2257760263159
be2_fci_6-31g_5.000_A.out:* ROOT 1 CI total energy = -29.2257572941079
be2_fci_6-31g_5.100_A.out:* ROOT 1 CI total energy = -29.2257367397319
be2_fci_6-31g_5.200_A.out:* ROOT 1 CI total energy = -29.2257152944106
be2_fci_6-31g_5.300_A.out:* ROOT 1 CI total energy = -29.2256936956757
be2_fci_6-31g_5.400_A.out:* ROOT 1 CI total energy = -29.2256725151707
be2_fci_6-31g_5.500_A.out:* ROOT 1 CI total energy = -29.2256521825575
be2_fci_6-31g_5.600_A.out:* ROOT 1 CI total energy = -29.2256330064081
be2_fci_6-31g_5.700_A.out:* ROOT 1 CI total energy = -29.2256151928252
be2_fci_6-31g_5.800_A.out:* ROOT 1 CI total energy = -29.2255988623269
be2_fci_6-31g_5.900_A.out:* ROOT 1 CI total energy = -29.2255840653183
be2_fci_6-31g_6.000_A.out:* ROOT 1 CI total energy = -29.2255707962738
be2_fci_6-31g_6.100_A.out:* ROOT 1 CI total energy = -29.2255590066328
be2_fci_6-31g_6.200_A.out:* ROOT 1 CI total energy = -29.2255486163380
be2_fci_6-31g_6.300_A.out:* ROOT 1 CI total energy = -29.2255395239277
be2_fci_6-31g_6.400_A.out:* ROOT 1 CI total energy = -29.2255316151255
be2_fci_6-31g_6.500_A.out:* ROOT 1 CI total energy = -29.2255247699107
be2_fci_6-31g_6.600_A.out:* ROOT 1 CI total energy = -29.2255188681184
be2_fci_6-31g_6.700_A.out:* ROOT 1 CI total energy = -29.2255137936650
be2_fci_6-31g_6.800_A.out:* ROOT 1 CI total energy = -29.2255094375408
be2_fci_6-31g_6.900_A.out:* ROOT 1 CI total energy = -29.2255056997400
be2_fci_6-31g_7.000_A.out:* ROOT 1 CI total energy = -29.2255024903089
be2_fci_6-31g_7.100_A.out:* ROOT 1 CI total energy = -29.2254997296976
be2_fci_6-31g_7.200_A.out:* ROOT 1 CI total energy = -29.2254973485845
be2_fci_6-31g_7.300_A.out:* ROOT 1 CI total energy = -29.2254952873300
be2_fci_6-31g_7.400_A.out:* ROOT 1 CI total energy = -29.2254934951891
be2_fci_6-31g_7.500_A.out:* ROOT 1 CI total energy = -29.2254919293907
be2_fci_6-31g_7.600_A.out:* ROOT 1 CI total energy = -29.2254905541657
be2_fci_6-31g_7.700_A.out:* ROOT 1 CI total energy = -29.2254893397894
be2_fci_6-31g_7.800_A.out:* ROOT 1 CI total energy = -29.2254882616727
be2_fci_6-31g_7.900_A.out:* ROOT 1 CI total energy = -29.2254872995360
be2_fci_6-31g_8.000_A.out:* ROOT 1 CI total energy = -29.2254864366749
be2_fci_6-31g_8.100_A.out:* ROOT 1 CI total energy = -29.2254856593233
be2_fci_6-31g_8.200_A.out:* ROOT 1 CI total energy = -29.2254849561130
be2_fci_6-31g_8.300_A.out:* ROOT 1 CI total energy = -29.2254843176196
be2_fci_6-31g_8.400_A.out:* ROOT 1 CI total energy = -29.2254837359900
be2_fci_6-31g_8.500_A.out:* ROOT 1 CI total energy = -29.2254832046383
be2_fci_6-31g_8.600_A.out:* ROOT 1 CI total energy = -29.2254827180001
be2_fci_6-31g_8.700_A.out:* ROOT 1 CI total energy = -29.2254822713372
be2_fci_6-31g_8.800_A.out:* ROOT 1 CI total energy = -29.2254818605812
be2_fci_6-31g_8.900_A.out:* ROOT 1 CI total energy = -29.2254814822109
be2_fci_6-31g_9.000_A.out:* ROOT 1 CI total energy = -29.2254811331546
be2_fci_6-31g_9.100_A.out:* ROOT 1 CI total energy = -29.2254808107136
be2_fci_6-31g_9.200_A.out:* ROOT 1 CI total energy = -29.2254805125017
be2_fci_6-31g_9.300_A.out:* ROOT 1 CI total energy = -29.2254802363964
be2_fci_6-31g_9.400_A.out:* ROOT 1 CI total energy = -29.2254799805006
be2_fci_6-31g_9.500_A.out:* ROOT 1 CI total energy = -29.2254797431115
be2_fci_6-31g_9.600_A.out:* ROOT 1 CI total energy = -29.2254795226950
be2_fci_6-31g_9.700_A.out:* ROOT 1 CI total energy = -29.2254793178647
be2_fci_6-31g_9.800_A.out:* ROOT 1 CI total energy = -29.2254791273648
be2_fci_6-31g_9.900_A.out:* ROOT 1 CI total energy = -29.2254789500548
be2_fci_6-31g_10.000_A.out:* ROOT 1 CI total energy = -29.2254787848972
be2_fci_6-31g_10.100_A.out:* ROOT 1 CI total energy = -29.2254786309469
be2_fci_6-31g_10.200_A.out:* ROOT 1 CI total energy = -29.2254784873412
be2_fci_6-31g_10.300_A.out:* ROOT 1 CI total energy = -29.2254783532923
be2_fci_6-31g_10.400_A.out:* ROOT 1 CI total energy = -29.2254782280791
be2_fci_6-31g_10.500_A.out:* ROOT 1 CI total energy = -29.2254781110419
be2_fci_6-31g_10.600_A.out:* ROOT 1 CI total energy = -29.2254780015758
be2_fci_6-31g_10.700_A.out:* ROOT 1 CI total energy = -29.2254778991259
be2_fci_6-31g_10.800_A.out:* ROOT 1 CI total energy = -29.2254778031827
be2_fci_6-31g_10.900_A.out:* ROOT 1 CI total energy = -29.2254777132779
be2_fci_6-31g_11.000_A.out:* ROOT 1 CI total energy = -29.2254776289811
be2_fci_6-31g_11.100_A.out:* ROOT 1 CI total energy = -29.2254775498957
be2_fci_6-31g_11.200_A.out:* ROOT 1 CI total energy = -29.2254774756568
be2_fci_6-31g_11.300_A.out:* ROOT 1 CI total energy = -29.2254774059277
be2_fci_6-31g_11.400_A.out:* ROOT 1 CI total energy = -29.2254773403980
be2_fci_6-31g_11.500_A.out:* ROOT 1 CI total energy = -29.2254772787811
be2_fci_6-31g_11.600_A.out:* ROOT 1 CI total energy = -29.2254772208120
be2_fci_6-31g_11.700_A.out:* ROOT 1 CI total energy = -29.2254771662461
be2_fci_6-31g_11.800_A.out:* ROOT 1 CI total energy = -29.2254771148569
be2_fci_6-31g_11.900_A.out:* ROOT 1 CI total energy = -29.2254770664346
be2_fci_6-31g_12.000_A.out:* ROOT 1 CI total energy = -29.2254770207851
be2_fci_6-31g_12.100_A.out:* ROOT 1 CI total energy = -29.2254769777283
be2_fci_6-31g_12.200_A.out:* ROOT 1 CI total energy = -29.2254769370973
be2_fci_6-31g_12.300_A.out:* ROOT 1 CI total energy = -29.2254768987371
be2_fci_6-31g_12.400_A.out:* ROOT 1 CI total energy = -29.2254768625038
be2_fci_6-31g_12.500_A.out:* ROOT 1 CI total energy = -29.2254768282636
be2_fci_6-31g_12.600_A.out:* ROOT 1 CI total energy = -29.2254767958921
be2_fci_6-31g_12.700_A.out:* ROOT 1 CI total energy = -29.2254767652739
be2_fci_6-31g_12.800_A.out:* ROOT 1 CI total energy = -29.2254767363010
be2_fci_6-31g_12.900_A.out:* ROOT 1 CI total energy = -29.2254767088731
be2_fci_6-31g_13.000_A.out:* ROOT 1 CI total energy = -29.2254766828969
be2_fci_6-31g_13.100_A.out:* ROOT 1 CI total energy = -29.2254766582852
be2_fci_6-31g_13.200_A.out:* ROOT 1 CI total energy = -29.2254766349566
be2_fci_6-31g_13.300_A.out:* ROOT 1 CI total energy = -29.2254766128353
be2_fci_6-31g_13.400_A.out:* ROOT 1 CI total energy = -29.2254765918503
be2_fci_6-31g_13.500_A.out:* ROOT 1 CI total energy = -29.2254765719354
be2_fci_6-31g_13.600_A.out:* ROOT 1 CI total energy = -29.2254765530287
be2_fci_6-31g_13.700_A.out:* ROOT 1 CI total energy = -29.2254765350723
be2_fci_6-31g_13.800_A.out:* ROOT 1 CI total energy = -29.2254765180119
be2_fci_6-31g_13.900_A.out:* ROOT 1 CI total energy = -29.2254765017969
be2_fci_6-31g_14.000_A.out:* ROOT 1 CI total energy = -29.2254764863796
be2_fci_6-31g_14.100_A.out:* ROOT 1 CI total energy = -29.2254764717156
be2_fci_6-31g_14.200_A.out:* ROOT 1 CI total energy = -29.2254764577629
be2_fci_6-31g_14.300_A.out:* ROOT 1 CI total energy = -29.2254764444826
be2_fci_6-31g_14.400_A.out:* ROOT 1 CI total energy = -29.2254764318378
be2_fci_6-31g_14.500_A.out:* ROOT 1 CI total energy = -29.2254764197939
be2_fci_6-31g_14.600_A.out:* ROOT 1 CI total energy = -29.2254764083185
be2_fci_6-31g_14.700_A.out:* ROOT 1 CI total energy = -29.2254763973813
be2_fci_6-31g_14.800_A.out:* ROOT 1 CI total energy = -29.2254763869534
be2_fci_6-31g_14.900_A.out:* ROOT 1 CI total energy = -29.2254763770081
be2_fci_6-31g_15.000_A.out:* ROOT 1 CI total energy = -29.2254763675198
be2_fci_6-31g_15.100_A.out:* ROOT 1 CI total energy = -29.2254763584649
be2_fci_6-31g_15.200_A.out:* ROOT 1 CI total energy = -29.2254763498208
be2_fci_6-31g_15.300_A.out:* ROOT 1 CI total energy = -29.2254763415662
be2_fci_6-31g_15.400_A.out:* ROOT 1 CI total energy = -29.2254763336813
be2_fci_6-31g_15.500_A.out:* ROOT 1 CI total energy = -29.2254763261475
be2_fci_6-31g_15.600_A.out:* ROOT 1 CI total energy = -29.2254763189467
be2_fci_6-31g_15.700_A.out:* ROOT 1 CI total energy = -29.2254763120623
be2_fci_6-31g_15.800_A.out:* ROOT 1 CI total energy = -29.2254763054786
be2_fci_6-31g_15.900_A.out:* ROOT 1 CI total energy = -29.2254762991806
be2_fci_6-31g_16.000_A.out:* ROOT 1 CI total energy = -29.2254762931541
be2_fci_6-31g_16.100_A.out:* ROOT 1 CI total energy = -29.2254762873860
be2_fci_6-31g_16.200_A.out:* ROOT 1 CI total energy = -29.2254762818637
be2_fci_6-31g_16.300_A.out:* ROOT 1 CI total energy = -29.2254762765752
be2_fci_6-31g_16.400_A.out:* ROOT 1 CI total energy = -29.2254762715093
be2_fci_6-31g_16.500_A.out:* ROOT 1 CI total energy = -29.2254762666554
be2_fci_6-31g_16.600_A.out:* ROOT 1 CI total energy = -29.2254762620033
be2_fci_6-31g_16.700_A.out:* ROOT 1 CI total energy = -29.2254762575435
be2_fci_6-31g_16.800_A.out:* ROOT 1 CI total energy = -29.2254762532671
be2_fci_6-31g_16.900_A.out:* ROOT 1 CI total energy = -29.2254762491654
be2_fci_6-31g_17.000_A.out:* ROOT 1 CI total energy = -29.2254762452303
be2_fci_6-31g_17.100_A.out:* ROOT 1 CI total energy = -29.2254762414541
be2_fci_6-31g_17.200_A.out:* ROOT 1 CI total energy = -29.2254762378296
be2_fci_6-31g_17.300_A.out:* ROOT 1 CI total energy = -29.2254762343498
be2_fci_6-31g_17.400_A.out:* ROOT 1 CI total energy = -29.2254762310082
be2_fci_6-31g_17.500_A.out:* ROOT 1 CI total energy = -29.2254762277985
be2_fci_6-31g_17.600_A.out:* ROOT 1 CI total energy = -29.2254762247147
be2_fci_6-31g_17.700_A.out:* ROOT 1 CI total energy = -29.2254762217513
be2_fci_6-31g_17.800_A.out:* ROOT 1 CI total energy = -29.2254762189030
be2_fci_6-31g_17.900_A.out:* ROOT 1 CI total energy = -29.2254762161647
be2_fci_6-31g_18.000_A.out:* ROOT 1 CI total energy = -29.2254762135314
be2_fci_6-31g_18.100_A.out:* ROOT 1 CI total energy = -29.2254762109987
be2_fci_6-31g_18.200_A.out:* ROOT 1 CI total energy = -29.2254762085622
be2_fci_6-31g_18.300_A.out:* ROOT 1 CI total energy = -29.2254762062178
be2_fci_6-31g_18.400_A.out:* ROOT 1 CI total energy = -29.2254762039614
be2_fci_6-31g_18.500_A.out:* ROOT 1 CI total energy = -29.2254762017893
be2_fci_6-31g_18.600_A.out:* ROOT 1 CI total energy = -29.2254761996980
be2_fci_6-31g_18.700_A.out:* ROOT 1 CI total energy = -29.2254761976840
be2_fci_6-31g_18.800_A.out:* ROOT 1 CI total energy = -29.2254761957439
be2_fci_6-31g_18.900_A.out:* ROOT 1 CI total energy = -29.2254761938750
be2_fci_6-31g_19.000_A.out:* ROOT 1 CI total energy = -29.2254761920740
be2_fci_6-31g_19.100_A.out:* ROOT 1 CI total energy = -29.2254761903383
be2_fci_6-31g_19.200_A.out:* ROOT 1 CI total energy = -29.2254761886651
be2_fci_6-31g_19.300_A.out:* ROOT 1 CI total energy = -29.2254761870518
be2_fci_6-31g_19.400_A.out:* ROOT 1 CI total energy = -29.2254761854960
be2_fci_6-31g_19.500_A.out:* ROOT 1 CI total energy = -29.2254761839955
be2_fci_6-31g_19.600_A.out:* ROOT 1 CI total energy = -29.2254761825479
be2_fci_6-31g_19.700_A.out:* ROOT 1 CI total energy = -29.2254761811511
be2_fci_6-31g_19.800_A.out:* ROOT 1 CI total energy = -29.2254761798031
be2_fci_6-31g_19.900_A.out:* ROOT 1 CI total energy = -29.2254761785021
be2_fci_6-31g_20.000_A.out:* ROOT 1 CI total energy = -29.2254761772460
""")



N = len(data)
start_doubled_index = 30	# The data starts at 3.0.  6.0 is the 30th element of the data
inc_trimer = []
for i in range(N):
	j = start_doubled_index + i*2
	if i<N and j<N:
		dist,  NNenergy = data[i]
		_,    NNNenergy = data[j]
		inc_trimer += [ (dist, 2*NNenergy+NNNenergy) ]
