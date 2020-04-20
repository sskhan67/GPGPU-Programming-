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
from   classes     import operator, composite_operator, commutator
from   function    import rule_distributor
from   parameters  import N, c1, e, f, d, p, r, q, s
from   parameters  import c2, k, l, h, a, b, i, j

virtual =  ['a', 'b', 'c', 'd']
occupied = ['i', 'j', 'k', 'l']


comutator = commutator(operator, composite_operator, rule_distributor)
for A in virtual:
	for I in occupied:
		for C in virtual:
			for K in occupied:
					for D in virtual:
						for L in occupied:
							for E in virtual:
								for M in occupied:
									for F in virtual:
										for N in occupied:
											g   = [+c1, [e, f, d], [i, j, a, b] ]
											t1 = [c2, k, [A, I] ]
											t2 = [c2, k, [C, K] ]
											t3 = [c2, k, [D, L] ]
											t4 = [c2, k, [E, M] ]
											t5 = [c2, k, [F, N] ]
											T1 = operator(t1)
											T2 = operator(t2)
											T3 = operator(t3)
											T4 = operator(t4)
											T5 = operator(t5)
											G  = operator( g )
											com1 = comutator(G,  T1 )
											for Z in range( len(com1) ):
												com2 = comutator(com1[Z], T2)
												for X in range( len(com2) ):
													com3 = comutator( com2[X], T3 )
													for Y in range( len(com3)  ):
														com4 = comutator( com3[Y], T4 )
														for W in range( len(com4) ):
															com5 = comutator(com4[W], T5 )
															if len(com5) != 0:
																print(com5[0].full_op)
