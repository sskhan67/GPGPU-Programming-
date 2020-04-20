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
from copy import deepcopy
from math import pi, sqrt, exp

"""\
General theory:

Each integral [overlap (S), kinetic energy (T), Coulomb (V)] over spherical (s) Gaussian orbitals or pair distributions
may be expressed in terms of a single-parameter kernel, generically K(U), whereby U alone carries all dependencies on
the nuclear centers.  Furthermore, U may be written such that it depends (nonlinearly) on three parameters which each
vary linearly with the 6 or 12 nuclear coordinates in question.

Since higher angular momenta orbitals are generated by differentiating functions with respect to these nuclear centers,
by the  chain rule, it suffices to know all of the full derivatives of K (up to finite order), and the partial derivatives
of U (with respect to the its three nonlinear parameters).  The target integrals are then a linear transformation of these
primitive quantities.

The big advantage of knowing this is that much of the math for handling higher angular momenta becomes independent
of the kernel K.

"""



max_L = 2	# the maximum angular momentum of any given basis function.

def resolve_derivs(orb):
	"""\
	The input and output are a representations of the same orbital, which may be written
		c * x^i * y^j * z^k * e^[-p * (r-P)^2]
		where r=<x,y,z>  and P=<Px,Py,Pz>
	It comes in as a nested tuple of the form:
		((Px,Py,Pz), p, c, (i,j,k))
	The output is a tuple of the form:
		((Px,Py,Pz), p, deriv_list)
	where deriv_list is of the form:
		[(c1, (i1,j1,k1)), (c2, (i2,j2,k2)), ... ]
	implicitly representing a sum, where each term is of the form
		c1 * (d/dPx)^i1 * (d/dPy)^j1 * (d/dPz)^k1 * e^[-p * (r-P)^2]
	and building an alternate representation identical to the input.
	"""
	P, p, c, el = orb
	c2 = c/(2*p)
	c4 = c/(4*(p**2))
	i,j,k = el
	L = i + j + k
	if L>max_L:  Exception("Anglular momentum exceeds hard-coded limit")
	if   L==0:
		return (P, p, [(c,el)] )
	elif L==1:
		return (P, p, [(c2, el)] )
	elif L==2:
		zero = (0,0,0)
		if 2 in (i, j, k):
			return (P, p, [(c4,el), (c2,zero)] )
		else:
			return (P, p, [(c4,el)] )
	else:
		raise Exception("Hard-coded angular momentum limit not aligned with actual implementation")



dx, dy, dz, dU = "dx", "dy", "dz", "dU"

def product_rule(deriv, terms):
	new = []
	for term in terms:
		for n,factor in enumerate(term):
			new_term = deepcopy(term)
			new_term[n] = factor + [deriv]
			new += [new_term]
	return new

def differentiate_product(derivs, n_factors):
	terms = []
	for _ in range(n_factors):  terms += [[]]
	for deriv in derivs:  terms = product_rule(deriv,terms)
	return terms

def chain_product_rule(deriv, terms):
	"""\
	terms is a list, where each element is a representation of a term in a sum.
	Each term has the form:
		[(d/dU)^n K(U)] * [(d/dx)U * (d/dz)(d/dy)U * ... ]
	where the product of derivatives of the variable U is relatively arbitrary, so
	long as the total number of such factors equals n, and the total number of
	differentiations of U must be the same in each term.
	Therefore, each element of terms holds a list of factors, where each factor is itself
	represented by a list of derivatives acted on a single U.  Summing the number of
	such factors gives the order n of the derivative of K.
	So the nesting is:
		list of terms
			each term is a list of U derivs
				each T deriv is a list of derivs to take on that U
	Calling this function updates this list as if the indicated derivative has
	been acted on the represented function.  It is presumed that U depends on three
	variables, x, y and z.
	"""
	new = []
	for term in terms:
		new += [ term + [[deriv]] ]
		for n,factor in enumerate(term):
			new_term = deepcopy(term)
			new_term[n] = factor + [deriv]
			new += [new_term]
	return new

def differentiate_kernel(derivs):
	"""\
	derivs is a list of arbitrary length containing only the values "dx", "dy", and "dz",
	representing derivatives with respect to x, y and z.  These are acted onto the function
	K(U), where U is taken to be an implicit function of x, y and z.  A list of terms of 
	the form described in the comments to chain_product_rule() is returned.
	"""
	terms = [[]]
	for deriv in derivs:  terms = chain_product_rule(deriv,terms)
	return terms



def mu(p, q):
	a = p + q
	return (p*q)/a

def R2(R):
	derivs = {}
	for i in range(3):
		for j in range(3):
			for k in range(3):
				derivs[i,j,k] = 0
	Rx, Ry, Rz = R
	derivs[0,0,0] = Rx**2 + Ry**2 + Rz**2
	derivs[1,0,0] = 2*Rx
	derivs[0,1,0] = 2*Ry
	derivs[0,0,1] = 2*Rz
	derivs[2,0,0] = 2
	derivs[0,2,0] = 2
	derivs[0,0,2] = 2
	return derivs

def U(p, P, q, Q):
	Px, Py, Pz = P
	Qx, Qy, Qz = Q
	Rx, Ry, Rz = Px-Qx, Py-Qy, Pz-Qz
	derivs = {}
	for key,value in R2(R).items()
		derivs[key] = mu(p,q) * value

def phase_derivs(derivs):
	phase = 1
	for n in range(len(derivs)):
		derivs[n] *= phase
		phase *= -1
	return derivs

def S_kernel(p, q, Uval):
	a = p + q
	derivs = []
	constant = sqrt((pi/a)**3) * exp(-Uval)		# this kernel is itself under differentiation with respect to -U
	for _ in range(max_L+1):
		derivs += [ constant ]
	return phase_derivs(derivs)			# changes derivatives with respect to -U into ones with respect to U

def T_raw(p, q, Uval):
	m = mu(p,q)
	derivs = [0 for _ in range(max_L+1)]
	derivs[0] = m * (3 - 2*Uval)
	derivs[1] = -2*m

def T_kernel(p, q, Uval):
	T_raw_derivs = T_raw(p, q, Uval)
	S_derivs = S_kernel(p, q, Uval)
	d = []
	derivs = []
	for _ in range(max_L+1):
		differentiate_product(d,2)
		result = 0
		for term in terms:
			oS = len(term[0])
			oT = len(term[1])
			result += S_derivs[oS] * T_raw_derivs[oT]
		derivs += [result]
		d += [dU]

def int_1e(gP, gQ, kernel):
	P, p, clP = resolve_derivs(gP)
	Q, q, clQ = resolve_derivs(gQ)
	U_derivs = U(p,P,q,Q)
	kernel_derivs = kernel(p, q, U[0,0,0])
	result = 0
	for cP,(iP,jP,kP) in clP:
		for cQ,(iQ,jQ,kQ) in clQ:
			derivs = []
			for _ in range(iP+iQ):  derivs += [dx]
			for _ in range(jP+jQ):  derivs += [dy]
			for _ in range(kP+kQ):  derivs += [dz]
			phase = (-1)**(iQ+jQ+kQ)
			terms = differentiate_kernel(derivs)
			for term in terms:
				subresult = kernel_derivs[len(term)]
				for Uderiv in term:
					if len(Uderiv)<3:
						i,j,k = 0,0,0
						for d in Uderiv:
							if d==dx:  i += 1
							if d==dy:  j += 1
							if d==dz:  k += 1
						subresult *= U_derivs[i,j,k]
				result += cP*cQ * subresult
	return result











def C(p, P, q, Q, derivs=[]):
	a = p + q
	Px,Py,Pz = P
	Qx,Qy,Qz = Q
	L = len(derivs)
	if   L==0:
		Cx = p*Px + q*Qx
		Cy = p*Py + q*Qy
		Cz = p*Pz + q*Qz
		return (Cx/a, Cy/a, Cz/a)
	elif L==1:
		d, = derivs
		if d==dPx:  return (p*Px/a, 0,      0     )
		if d==dPy:  return (0,      p*Py/a, 0     )
		if d==dPz:  return (0,      0,      p*Pz/a)
		if d==dQx:  return (q*Qx/a, 0,      0     )
		if d==dQy:  return (0,      q*Qy/a, 0     )
		if d==dQz:  return (0,      0,      q*Qz/a)
	else:
		return (0,0,0)