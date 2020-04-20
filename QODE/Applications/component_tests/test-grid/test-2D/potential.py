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

from parameters  import *

def HO_centered_at(x0):		# returns a shifted potential function
    kF = 1.			# force constant
    def v(x):
        x = x - x0
        return kF * x**2 / 2.
    return v


def HO_2D_centered_at(x0,y0):         # returns a shifted potential function
    kF1 = 1.                        # force constant
    kF2 = 1.
    def v(x,y):
        x = x - x0
        y = y - y0
        return ( (kF1 * x**2 / 2.) + (kF2 * y**2 / 2.) )
    return v


def two_particle_barrier_HO(height,radius,x0,kF): # returns a potential function
    def v(x1,x2):
        x2  = x2 - x0
        HO =  kF * x2**2 / 2.
        if abs(x1 - x2) < radius:
           interaction = height           
        else:
           interaction = 0
        return interaction + HO
    return v

def barrier_2D_at(x0, y0):
        L = 3
        H = 20
        end_barrier = x0 + L/2
        beg_barrier = x0 - L/2
        def v(x, y):
                if   x <=  beg_barrier:
                        return 0
                elif x >=  end_barrier:
                        return 0
                else:
                        return H
        return v

def barrier_at(x0):
	L = 3
	H = 20
	end_barrier = x0 + L/2
	beg_barrier = x0 - L/2
	def v(x):
		if   x <=  beg_barrier:
			return 0
		elif x >=  end_barrier:
			return 0
		else:
			return H
	return v 

def wall_at(x0):
	beginning = 3 + x0
	Height    = 5
	def v(x):
		if x >= beginning:
			return 0
		else:
			return Height
	return v	

def barrier_moved_from(x0):
	d = x0 + 5
	L = 3
	H = 5
	end_barrier = d + L/2
	beg_barrier = d - L/2
	def v(x):
		if   x <=  beg_barrier:
			return 0
		elif x >=  end_barrier:
			return 0
		else:
			return H
	return v 
