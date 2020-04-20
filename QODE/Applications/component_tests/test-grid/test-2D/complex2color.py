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
import  math
import cmath


color_points = [
  (0,       (0.18000,0.01800,0.01800)) \
, (0.16667, (0.17000,0.08500,0.00000)) \
, (0.33333, (0.12500,0.12500,0.00000)) \
, (0.50000, (0.04000,0.12000,0.04000)) \
, (0.58333, (0.04500,0.12000,0.08000)) \
, (0.66667, (0.06429,0.09000,0.15429)) \
, (0.75000, (0.10500,0.06000,0.16500)) \
, (0.83333, (0.13750,0.01250,0.13750)) \
, (1,       (0.18000,0.01800,0.01800)) \
]


def complex2color_log(log_min=-4,log_max=0):
	Delta = log_max - log_min
	def map(v):
		r,phi = cmath.polar(v)
		phi /= 2*math.pi
		if phi<0:  phi += 1
		#
		if r>0:
			r = math.log10(r)
			r = (r-log_min)/Delta
			if r>1:  r = 1
			if r<0:  r = 0
		else:
			r = 0
		r *= 5
		#
		i = 0
		while phi>color_points[i][0]:  i += 1
		mix_high = (phi - color_points[i-1][0]) / (color_points[i][0] - color_points[i-1][0])
		mix_low  = 1 - mix_high
		Hr,Hg,Hb  = color_points[i][1]
		Lr,Lg,Lb  = color_points[i-1][1]
		R = mix_high*Hr + mix_low*Lr
		G = mix_high*Hg + mix_low*Lg
		B = mix_high*Hb + mix_low*Lb
		#
		R = (math.floor(255*r*R))
		G = (math.floor(255*r*G))
		B = (math.floor(255*r*B))
		return (255,R,G,B)
	return map

def complex2color(max_val):
	def map(v):
		r,phi = cmath.polar(v)
		phi /= 2*math.pi
		if phi<0:  phi += 1
		#
		r = r / max_val
		if r>1:  r = 1
		if r<0:  r = 0
		r *= 5
		#
		i = 0
		while phi>color_points[i][0]:  i += 1
		mix_high = (phi - color_points[i-1][0]) / (color_points[i][0] - color_points[i-1][0])
		mix_low  = 1 - mix_high
		Hr,Hg,Hb  = color_points[i][1]
		Lr,Lg,Lb  = color_points[i-1][1]
		R = mix_high*Hr + mix_low*Lr
		G = mix_high*Hg + mix_low*Lg
		B = mix_high*Hb + mix_low*Lb
		#
		R = int(math.floor(255*r*R))
		G = int(math.floor(255*r*G))
		B = int(math.floor(255*r*B))
		return (255,R,G,B)
	return map
