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
from letters               import *
from equations_transformer import normal_CC_to_TTCC 
from code_generation       import code
from latex_syntax          import latex

'''
Things to keep in mind about the class input:

 * This function only accepts its inputs as a list of lists. Although, the first and last arguments are not contained in a list.
 * The first argument is the prefactor, so it has to be a number.
 * The last argument is one of two strings('unrestricted' or 'restricted'). It controls the restrictions on the spacing limits: 'unrestricted' == -oo to +oo; 'restricted' == restrictions given by get_restrictions().
 * No need to input the summation indices of orbitals. The code will deduce them automatically.
 * The indices of a single integral or amplitude are contained in its own list.

Example of input
for the term in conventional coupled cluster(latex syntax):

\frac{1}{4}\sum_{klcd}V^{kl}_{cd}t^{c}_{i}t^{a}_{k}t^{d}_{n}t^{b}_{l}

the input syntax is:

normal_term = [1/4, [k, l, c, d], [c, i], [a, k], [d, j], [b, l], 'unrestricted' ]

'''


# Example inputs

example1 = [1/2, [k, a, c, d], [c, d, k, i], 'unrestricted' ]
example2 = [-1/2, [k, l, c, d], [c, d, k, i], [a, l], 'unrestricted' ]
example3 = [1, [a, c], [c, i], 'restricted']



##########################
# energy equation
E1 = [+1, [i, a], [a,i], 'restricted']
E2 = [+1/4, [i, j, a, b], [a, b, i, j], 'restricted']
E3 = [+1/2, [i, j, a, b], [a, i], [b, j], 'restricted']

# T1 equation
S1 = [+1, [a, i], 'restricted']
S2 = [+1, [a, c], [c, i], 'restricted']
S3 = [-1, [k, i], [a, k], 'restricted']
S4 = [+1, [k, a, c, i], [c, k], 'restricted']
S5 = [+1, [k, c], [a, c, i, k], 'restricted']
S6 = [+1/2, [k, a, c, d], [c, d, k, i], 'restricted']


S7  = [-1/2, [k, l, c, i], [c, a, k, l], 'restricted']
S8  = [-1, [k, c], [c, i], [a, k], 'restricted']
S9  = [-1, [k, l, c, i], [c, k], [a, l], 'restricted']
S10 = [-1, [k, a, c, d], [c, k], [d, i], 'restricted']


S11 = [-1, [k, l, c, d], [c, k], [d, i], [a, l], 'restricted']
S12 = [+1, [k, l, c, d], [c, k], [d, a, l, i], 'restricted']
S13 = [-1/2, [k, l, c, d], [c, d, k, i], [a, l], 'restricted']
S14 = [-1/2, [k, l, c, d], [c, a, k, l], [d, i], 'restricted']
##########################




###########################
# T2 equation
D1 = [1, [a, b, i, j], 'restricted']
D2 = [1, [b, c], [a, c, i, j], 'restricted']
D3 = [-1, [a, c], [b, c, i, j], 'restricted']
D4 = [-1, [k, j], [a, b, i, k], 'restricted']
D5 = [+1, [k, i], [a, b, j, k], 'restricted']
D6 = [1/2, [k, l, i, j], [a, b, k, l], 'restricted']


D7  = [1/2, [a, b, c, d], [c, d, i, j], 'restricted']
D8  = [ 1, [k, b, c, j], [a, c, i, k], 'restricted']
D9  = [-1, [k, b, c, i], [a, c, j, k], 'restricted']
D10 = [-1, [k, a, c, j], [b, c, i, k], 'restricted']
D11 = [ 1, [k, a, c, i], [b, c, j, k], 'restricted']

D12 = [1, [a, b, c, j], [c, i], 'restricted']
D13 = [-1, [a, b, c, i], [c, j], 'restricted']  

D14 = [-1, [k, b, i, j],[a, k], 'restricted']
D15 = [1, [k, a, i, j],[b, k], 'restricted']

D16  = [+1/2, [k, l, c, d], [a, c, i, k], [d, b, l, j], 'restricted']
D17  = [-1/2, [k, l, c, d], [b, c, i, k], [d, a, l, j], 'restricted']
D18  = [-1/2, [k, l, c, d], [a, c, j, k], [d, b, l, i], 'restricted']
D19  = [+1/2, [k, l, c, d], [b, c, j, k], [d, a, l, i], 'restricted']

D20 = [1/4, [k, l, c, d], [c, d, i, j], [a, b, k, l], 'restricted']

D21 = [-1/2, [k, l, c, d], [a, c, i, j], [b, d, k ,l], 'restricted']
D22 = [+1/2, [k, l, c, d], [b, c, i, j], [a, d, k ,l], 'restricted']

D23 = [-1/2, [k, l, c, d], [a, b, i, k], [c, d, j, l], 'restricted']
D24 = [+1/2, [k, l, c, d], [a, b, j, k], [c, d, i, l], 'restricted']

D25 = [+1/2, [k, l, i, j], [a, k], [b, l], 'restricted']
D26 = [-1/2, [k, l, i, j], [b, k], [a, l], 'restricted']


D27 = [+1/2, [a, b, c, d], [c, i], [d, j], 'restricted']
D28 = [-1/2, [a, b, c, d], [c, j], [d, i], 'restricted']

D29   = [-1, [k, b, i, c], [a, k], [c, j], 'restricted']
D30   = [+1, [k, a, i, c], [b, k], [c, j], 'restricted']
D31   = [+1, [k, b, j, c], [a, k], [c, i], 'restricted']
D32   = [-1, [k, a, j, c], [b, k], [c, i], 'restricted']

D33 = [ 1, [k, c], [a, k], [b, c, i, j], 'restricted']
D34 = [-1, [k, c], [b, k], [a, c, i, j], 'restricted']


D35   = [ 1, [k, c], [c, i], [a, b, j, k], 'restricted']
D36 = [-1, [k, c], [c, j], [a, b, i, k], 'restricted']

D37 = [-1, [k, l, c, i], [c, k], [a, b, l, j], 'restricted']
D38 = [+1, [k, l, c, j], [c, k], [a, b, l, i], 'restricted']

D39 = [ 1, [k, a, c, d], [c, k], [d, b, i, j], 'restricted']
D40 = [-1, [k, b, c, d], [c, k], [d, a, i, j], 'restricted']

D41 = [+1, [a, k, d, c], [d, i], [b, c, j, k], 'restricted']
D42 = [-1, [b, k, d, c], [d, i], [a, c, j, k], 'restricted']
D43 = [-1, [a, k, d, c], [d, j], [b, c, i, k], 'restricted']
D44 = [+1, [b, k, d, c], [d, j], [a, c, i, k], 'restricted']

D45 = [1,  [k, l, i, c], [a, l], [b, c, j, k], 'restricted']
D46 = [-1, [k, l, i, c], [b, l], [a, c, j, k], 'restricted']
D47 = [-1, [k, l, j, c], [a, l], [b, c, i, k], 'restricted']
D48 = [1,  [k, l, j, c], [b, l], [a, c, i, k], 'restricted']


D49 = [ 1/2, [k, l, c, j], [c, i], [a, b, k, l], 'restricted']
D50 = [-1/2, [k, l, c, i], [c, j], [a, b, k, l], 'restricted']


D51 = [-1/2, [k, b, c, d], [a, k], [c, d, i, j], 'restricted']
D52 = [+1/2, [k, a, c, d], [b, k], [c, d, i, j], 'restricted']


D53 = [-1/2, [k, b, c, d], [c, i], [a, k], [d, j], 'restricted']
D54 = [+1/2, [k, a, c, d], [c, i], [b, k], [d, j], 'restricted']
D55 = [+1/2, [k, b, c, d], [c, j], [a, k], [d, i], 'restricted']
D56 = [-1/2, [k, a, c, d], [c, j], [b, k], [d, i], 'restricted']

D57 = [ 1/2, [k, l, c, j], [c, i], [a, k], [b, l], 'restricted']
D58 = [-1/2, [k, l, c, j], [c, i], [b, k], [a, l], 'restricted']
D59 = [-1/2, [k, l, c, i], [c, j], [a, k], [b, l], 'restricted']
D60 = [ 1/2, [k, l, c, i], [c, j], [b, k], [a, l], 'restricted']


D61 = [-1, [k, l, c, d], [c, k], [d, i], [a, b, l, j], 'restricted']
D62 = [+1, [k, l, c, d], [c, k], [d, j], [a, b, l, i], 'restricted']


D63 = [-1, [k, l, c, d], [c, k], [a, l], [d, b, i, j], 'restricted']
D64 = [+1, [k, l, c, d], [c, k], [b, l], [d, a, i, j], 'restricted']


D65  = [ 1/4, [k, l, c, d], [c, i], [d, j], [a, b, k, l], 'restricted']
D66  = [-1/4, [k, l, c, d], [c, j], [d, i], [a, b, k, l], 'restricted']


D67 = [ 1/4, [k, l, c, d], [a, k], [b, l], [c, d, i, j], 'restricted']
D68 = [-1/4, [k, l, c, d], [b, k], [a, l], [c, d, i, j], 'restricted']


D69 = [ 1,   [k, l, c, d], [c, i], [b, l], [a, d, k, j], 'restricted']
D70 = [-1,   [k, l, c, d], [c, i], [a, l], [b, d, k, j], 'restricted']
D71 = [-1,   [k, l, c, d], [c, j], [b, l], [a, d, k, i], 'restricted']
D72 = [ 1,   [k, l, c, d], [c, j], [a, l], [b, d, k, i], 'restricted']


D73 = [ 1/4, [k, l, c, d], [c, i], [a, k], [d, j], [b, l], 'restricted']
D74 = [-1/4, [k, l, c, d], [c, i], [b, k], [d, j], [a, l], 'restricted']
D75 = [-1/4, [k, l, c, d], [c, j], [a, k], [d, i], [b, l], 'restricted']
D76 = [ 1/4, [k, l, c, d], [c, j], [b, k], [d, i], [a, l], 'restricted']
#########################################################################

Example_equation= [example1, example2, example3]
Energy_equation = [E1, E2, E3]
T1_equation     = [S1,  S2,  S3,  S4,  S5,  S6,  S7,  S8,  S9,  S10, S11, S12, S13, S14]
T2_equation     = [D1,  D2,  D3,  D4,  D5,  D6,  D7,  D8,  D9,  D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20,
                   D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40,
                   D41, D42, D43, D44, D45, D46, D47, D48, D49, D50, D51, D52, D53, D54, D55, D56, D57, D58, D59, D60,
                   D61, D62, D63, D64, D65, D66, D67, D68, D69, D70, D71, D72, D73, D74, D75, D76]

file_name = "generated_CC_equation_code.py"
code(T2_equation, file_name, 'translationally transformed CC') 
#latex(T2_equation)






