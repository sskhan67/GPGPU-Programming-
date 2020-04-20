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
	dist   = float(".".join(fields[0].split("_")[1].split(".")[:2]))
	energy = float(fields[-1])
	return dist, energy-3*atom.energy

data = X(
parse,
"""\
Be3_3.8.out:      * CCSD(T) total energy                  =  -43.837655710567866
Be3_3.9.out:      * CCSD(T) total energy                  =  -43.837989016717970
Be3_4.0.out:      * CCSD(T) total energy                  =  -43.838249944564488
Be3_4.1.out:      * CCSD(T) total energy                  =  -43.838448727848821
Be3_4.2.out:      * CCSD(T) total energy                  =  -43.838595046020053
Be3_4.3.out:      * CCSD(T) total energy                  =  -43.838697756145777
Be3_4.4.out:      * CCSD(T) total energy                  =  -43.838764765637201
Be3_4.5.out:      * CCSD(T) total energy                  =  -43.838803001184004
Be3_4.6.out:      * CCSD(T) total energy                  =  -43.838818439455380
Be3_4.7.out:      * CCSD(T) total energy                  =  -43.838816174652841
Be3_4.8.out:      * CCSD(T) total energy                  =  -43.838800504084475
Be3_4.9.out:      * CCSD(T) total energy                  =  -43.838775019076770
Be3_5.0.out:      * CCSD(T) total energy                  =  -43.838742693633698
Be3_5.1.out:      * CCSD(T) total energy                  =  -43.838705965801445
Be3_5.2.out:      * CCSD(T) total energy                  =  -43.838666810114816
Be3_5.3.out:      * CCSD(T) total energy                  =  -43.838626800884953
Be3_5.4.out:      * CCSD(T) total energy                  =  -43.838587167399105
Be3_5.5.out:      * CCSD(T) total energy                  =  -43.838548842105880
Be3_5.6.out:      * CCSD(T) total energy                  =  -43.838512503648936
Be3_5.7.out:      * CCSD(T) total energy                  =  -43.838478615388560
Be3_5.8.out:      * CCSD(T) total energy                  =  -43.838447460924371
Be3_5.9.out:      * CCSD(T) total energy                  =  -43.838419176317871
Be3_6.0.out:      * CCSD(T) total energy                  =  -43.838393780006633
Be3_6.1.out:      * CCSD(T) total energy                  =  -43.838371198904049
Be3_6.2.out:      * CCSD(T) total energy                  =  -43.838351293891122
Be3_6.3.out:      * CCSD(T) total energy                  =  -43.838333879206473
Be3_6.4.out:      * CCSD(T) total energy                  =  -43.838318741018533
Be3_6.5.out:      * CCSD(T) total energy                  =  -43.838305651518461
Be3_6.6.out:      * CCSD(T) total energy                  =  -43.838294381156622
Be3_6.7.out:      * CCSD(T) total energy                  =  -43.838284706612690
Be3_6.8.out:      * CCSD(T) total energy                  =  -43.838276417577198
Be3_6.9.out:      * CCSD(T) total energy                  =  -43.838269320578817
Be3_7.0.out:      * CCSD(T) total energy                  =  -43.838263241283173
Be3_7.1.out:      * CCSD(T) total energy                  =  -43.838258025428416
Be3_7.2.out:      * CCSD(T) total energy                  =  -43.838253538415536
Be3_7.3.out:      * CCSD(T) total energy                  =  -43.838249664502619
Be3_7.4.out:      * CCSD(T) total energy                  =  -43.838246305271809
Be3_7.5.out:      * CCSD(T) total energy                  =  -43.838243377800637
Be3_7.6.out:      * CCSD(T) total energy                  =  -43.838240812831970
Be3_7.7.out:      * CCSD(T) total energy                  =  -43.838238552930008
Be3_7.8.out:      * CCSD(T) total energy                  =  -43.838236550636431
Be3_7.9.out:      * CCSD(T) total energy                  =  -43.838234766914553
Be3_8.0.out:      * CCSD(T) total energy                  =  -43.838233169697261
Be3_8.1.out:      * CCSD(T) total energy                  =  -43.838231732673002
Be3_8.2.out:      * CCSD(T) total energy                  =  -43.838230434122721
Be3_8.3.out:      * CCSD(T) total energy                  =  -43.838229256137730
Be3_8.4.out:      * CCSD(T) total energy                  =  -43.838228183855357
Be3_8.5.out:      * CCSD(T) total energy                  =  -43.838227204833750
Be3_8.6.out:      * CCSD(T) total energy                  =  -43.838226308611787
Be3_8.7.out:      * CCSD(T) total energy                  =  -43.838225486312197
Be3_8.8.out:      * CCSD(T) total energy                  =  -43.838224730329031
Be3_8.9.out:      * CCSD(T) total energy                  =  -43.838224034104989
Be3_9.0.out:      * CCSD(T) total energy                  =  -43.838223391928850
Be3_9.1.out:      * CCSD(T) total energy                  =  -43.838222798797048
Be3_9.2.out:      * CCSD(T) total energy                  =  -43.838222250270952
Be3_9.3.out:      * CCSD(T) total energy                  =  -43.838221742485786
Be3_9.4.out:      * CCSD(T) total energy                  =  -43.838221271872548
Be3_9.5.out:      * CCSD(T) total energy                  =  -43.838220835326986
Be3_9.6.out:      * CCSD(T) total energy                  =  -43.838220430008697
Be3_9.7.out:      * CCSD(T) total energy                  =  -43.838220053364246
Be3_9.8.out:      * CCSD(T) total energy                  =  -43.838219703080718
Be3_9.9.out:      * CCSD(T) total energy                  =  -43.838219377059040
""")
