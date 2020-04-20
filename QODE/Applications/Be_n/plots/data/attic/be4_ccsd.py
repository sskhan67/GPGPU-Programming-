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
	return dist, energy-4*atom.energy

data = X(
parse,
"""\
Be4_3.8.out:      * CCSD total energy                     =  -58.449494949969178
Be4_3.9.out:      * CCSD total energy                     =  -58.450090564500222
Be4_4.0.out:      * CCSD total energy                     =  -58.450563529548610
Be4_4.1.out:      * CCSD total energy                     =  -58.450931171057412
Be4_4.2.out:      * CCSD total energy                     =  -58.451209681318979
Be4_4.3.out:      * CCSD total energy                     =  -58.451413788933309
Be4_4.4.out:      * CCSD total energy                     =  -58.451556612512292
Be4_4.5.out:      * CCSD total energy                     =  -58.451649642726636
Be4_4.6.out:      * CCSD total energy                     =  -58.451702808038576
Be4_4.7.out:      * CCSD total energy                     =  -58.451724590832860
Be4_4.8.out:      * CCSD total energy                     =  -58.451722168399968
Be4_4.9.out:      * CCSD total energy                     =  -58.451701560932150
Be4_5.0.out:      * CCSD total energy                     =  -58.451667774820791
Be4_5.1.out:      * CCSD total energy                     =  -58.451624934672154
Be4_5.2.out:      * CCSD total energy                     =  -58.451576401125649
Be4_5.3.out:      * CCSD total energy                     =  -58.451524874385811
Be4_5.4.out:      * CCSD total energy                     =  -58.451472484661281
Be4_5.5.out:      * CCSD total energy                     =  -58.451420871765627
Be4_5.6.out:      * CCSD total energy                     =  -58.451371255834921
Be4_5.7.out:      * CCSD total energy                     =  -58.451324501029831
Be4_5.8.out:      * CCSD total energy                     =  -58.451281173878584
Be4_5.9.out:      * CCSD total energy                     =  -58.451241594486874
Be4_6.0.out:      * CCSD total energy                     =  -58.451205888262606
Be4_6.1.out:      * CCSD total energy                     =  -58.451174026148806
Be4_6.2.out:      * CCSD total energy                     =  -58.451145864899992
Be4_6.3.out:      * CCSD total energy                     =  -58.451121180541264
Be4_6.4.out:      * CCSD total energy                     =  -58.451099696890367
Be4_6.5.out:      * CCSD total energy                     =  -58.451081109264841
Be4_6.6.out:      * CCSD total energy                     =  -58.451065102850663
Be4_6.7.out:      * CCSD total energy                     =  -58.451051367665428
Be4_6.8.out:      * CCSD total energy                     =  -58.451039608246937
Be4_6.9.out:      * CCSD total energy                     =  -58.451029550744444
Be4_7.0.out:      * CCSD total energy                     =  -58.451020947268475
Be4_7.1.out:      * CCSD total energy                     =  -58.451013577398029
Be4_7.2.out:      * CCSD total energy                     =  -58.451007248374047
Be4_7.3.out:      * CCSD total energy                     =  -58.451001794139067
Be4_7.4.out:      * CCSD total energy                     =  -58.450997073245020
Be4_7.5.out:      * CCSD total energy                     =  -58.450992966504216
Be4_7.6.out:      * CCSD total energy                     =  -58.450989374388840
Be4_7.7.out:      * CCSD total energy                     =  -58.450986214356881
Be4_7.8.out:      * CCSD total energy                     =  -58.450983418353076
Be4_7.9.out:      * CCSD total energy                     =  -58.450980930442270
Be4_8.0.out:      * CCSD total energy                     =  -58.450978704789577
Be4_8.1.out:      * CCSD total energy                     =  -58.450976703809879
Be4_8.2.out:      * CCSD total energy                     =  -58.450974896683896
Be4_8.3.out:      * CCSD total energy                     =  -58.450973257968990
Be4_8.4.out:      * CCSD total energy                     =  -58.450971766643178
Be4_8.5.out:      * CCSD total energy                     =  -58.450970405185188
Be4_8.6.out:      * CCSD total energy                     =  -58.450969158888164
Be4_8.7.out:      * CCSD total energy                     =  -58.450968015312121
Be4_8.8.out:      * CCSD total energy                     =  -58.450966963829806
Be4_8.9.out:      * CCSD total energy                     =  -58.450965995293679
Be4_9.0.out:      * CCSD total energy                     =  -58.450965101760360
Be4_9.1.out:      * CCSD total energy                     =  -58.450964276275847
Be4_9.2.out:      * CCSD total energy                     =  -58.450963512709336
Be4_9.3.out:      * CCSD total energy                     =  -58.450962805621884
Be4_9.4.out:      * CCSD total energy                     =  -58.450962150160493
Be4_9.5.out:      * CCSD total energy                     =  -58.450961541975147
Be4_9.6.out:      * CCSD total energy                     =  -58.450960977148640
Be4_9.7.out:      * CCSD total energy                     =  -58.450960452146042
Be4_9.8.out:      * CCSD total energy                     =  -58.450959963764731
Be4_9.9.out:      * CCSD total energy                     =  -58.450959509096670
""")
