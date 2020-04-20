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
import sys

if __name__ == '__main__':
	print('input file =', sys.argv[1])
	f = open(sys.argv[1], 'r')
	lines = f.readlines()
	f.close()

	mol_buf = []
	e_buf = ['###']
	t_buf = ['###']
	n_buf = ['###']

	for line in lines:
		if 'mol' in line:
			print(line)
			mol_buf += [line.strip('+ \n')] 
		elif 'TOTAL ENERGY' in line:
			print(line)
			mol_buf += [line.split()[-1]]
		elif 'ANALYTICAL SOLUTION' in line:
			print(line)
			e_buf += ['###']
			t_buf += ['###']
			n_buf += ['','###']
			mol_buf += [line.split()[-1]] 
		elif 'E_CCSD' in line and 'Projection' not in line:
			print(line.strip())
			e_buf += [ line.split()[-1] ]
		elif 'Relative step norm' in line:
			print(line.strip())
			n_buf += [ line.split()[-1] ]
		elif 'Cycle Cumulative Time' in line:
			print(line.strip())
			t_buf += [ line.split()[-1] ]

	e_buf = e_buf[:-1]
	t_buf = t_buf[:-1]
	n_buf = n_buf[:-1]

	print('len e_buf =', len(e_buf))
	print('len t_buf =', len(t_buf))
	print('len n_buf =', len(n_buf))

	output = []
	for i in range(len(e_buf)):
		output += [[ t_buf[i] , e_buf[i], n_buf[i] ]]

	#for item in output:
	#	print(item)

	ct = 0
	n_item = 0
	for item in output:
		if item == ['###','###','###']:
			cycle = 1
			output[n_item] = ['###' , mol_buf[ct*3] , mol_buf[ct*3+2] , mol_buf[ct*3+1] ] 
			ct += 1
		else:
			output[n_item] = [str(cycle)] + item
			cycle += 1
		n_item += 1				

	#for item in output:
	#	print(item)

	output_buf = 'Cycle	ElapsedTime(sec)	E_CCSD(a.u.)	RelNorm\n'
	for item in output:
		output_buf += '\t'.join(item) + '\n'
	print(output_buf)
	f=open(sys.argv[1].split('.')[0] + '.dat', 'w')
	f.write(output_buf)
	f.close()



