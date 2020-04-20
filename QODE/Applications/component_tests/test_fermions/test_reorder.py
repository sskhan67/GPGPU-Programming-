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
from qode.fermion_field         import occ_strings
from qode.fermion_field.state   import configuration, state, dot, op_string, create, annihilate
from qode.fermion_field.reorder import *



all_occ_strings = occ_strings.all_occ_strings(8,4)

occ = [0,2,4,6]
vrt = [1,3,5,7]





print("##################   p* q   ##################")
error = 0.

reorder = reorder_CpAq(vrt)
for p in occ+vrt:
	for q in occ+vrt:
		reordered = reorder(p,q)
		#
		print("\n{}* {}  \n   =".format(p,q))
		printable = ""
		sign_on = False
		for term in reordered:
			phase = term[-1]
			if sign_on:
				if phase>0:  printable += "   +   "
				else:        printable += "   -   "
				if len(term)==1:  printable += ("{:+d}".format(phase))[1:]
			else:
				if len(term)==1:  printable +=  "{: d}".format(phase)
				else:
					if phase>0:  printable += "  "
					else:        printable += "- "
			sign_on = True
			for op in term[:-1]:  printable += str(op)
		print(printable+"\n")
		for i in range(len(all_occ_strings)):
			printable = str(configuration(all_occ_strings[i]))
			psi_in = state(all_occ_strings,i)
			#
			control = op_string(create(p),annihilate(q))(psi_in)
			experiment = state(all_occ_strings)
			for term in reordered:
				phase = term[-1]
				temp = psi_in
				for t in reversed(term[:-1]):  temp = t(temp)
				experiment.increment(temp,phase)
			#
			control.increment(experiment,-1)
			err = dot(control,control)
			printable += "         error = {}".format(err)
			print(printable)
			error += err
			if err>1e-10:  raise Exception("inequivalent operators")

print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))



print("##################   p* q* r  s    ##################")
error = 0.

reorder = reorder_CpCqArAs(vrt)
for p in occ+vrt:
	for q in occ+vrt:
		for r in occ+vrt:
			for s in occ+vrt:
				reordered = reorder(p,q,r,s)
				#
				print("\n{}* {}* {}  {}  \n   =".format(p,q,r,s))
				printable = ""
				sign_on = False
				for term in reordered:
					phase = term[-1]
					if sign_on:
						if phase>0:  printable += "   +   "
						else:        printable += "   -   "
						if len(term)==1:  printable += ("{:+d}".format(phase))[1:]
					else:
						if len(term)==1:  printable +=  "{: d}".format(phase)
						else:
							if phase>0:  printable += "  "
							else:        printable += "- "
					sign_on = True
					for op in term[:-1]:  printable += str(op)
				print(printable+"\n")
				for i in range(len(all_occ_strings)):
					printable = str(configuration(all_occ_strings[i]))
					psi_in = state(all_occ_strings,i)
					#
					control = op_string(create(p),create(q),annihilate(r),annihilate(s))(psi_in)
					experiment = state(all_occ_strings)
					for term in reordered:
						phase = term[-1]
						temp = psi_in
						for t in reversed(term[:-1]):  temp = t(temp)
						experiment.increment(temp,phase)
					#
					control.increment(experiment,-1)
					err = dot(control,control)
					printable += "         error = {}".format(err)
					print(printable)
					error += err
					if err>1e-10:  raise Exception("inequivalent operators")

print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))





all_e = []
for i in occ:
	for a in vrt:
		all_e += [excitation(a,i)]
all_d = []
for i in occ:
	for a in vrt:
		all_d += [deexcitation(a,i)]
all_fo = []
for i in occ:
	for j in occ:
		all_fo += [occ_rearrange(i,j)]
all_fv = []
for a in vrt:
	for b in vrt:
		all_fv += [vrt_rearrange(a,b)]



def test_trans_reordering(trans_list1,trans_list2,all_occ_strings):
	error = 0.
	for t1 in trans_list1:
		for t2 in trans_list2:
			reordered = reorder_t_t(t1,t2)
			#
			print("\n" + str(t1)+str(t2) + "\n   =")
			printable = ""
			sign_on = False
			for term in reordered:
				phase = term[-1]
				if sign_on:
					if phase>0:  printable += "   +   "
					else:        printable += "   -   "
					if len(term)==1:  printable += ("{:+d}".format(phase))[1:]
				else:
					if len(term)==1:  printable +=  "{: d}".format(phase)
					else:
						if phase>0:  printable += "  "
						else:        printable += "- "
				sign_on = True
				for op in term[:-1]:  printable += str(op)
			print(printable+"\n")
			for i in range(len(all_occ_strings)):
				printable = str(configuration(all_occ_strings[i]))
				psi_in = state(all_occ_strings,i)
				#
				control = t1(t2(psi_in))
				experiment = state(all_occ_strings)
				for term in reordered:
					phase = term[-1]
					temp = psi_in
					for t in reversed(term[:-1]):  temp = t(temp)
					experiment.increment(temp,phase)
				#
				control.increment(experiment,-1)
				err = dot(control,control)
				printable += "         error = {}".format(err)
				print(printable)
				error += err
				if err>1e-10:  raise Exception("inequivalent operators")
	return error

print("##################   E  E    ##################")
error = test_trans_reordering(all_e,all_e,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   E  D    ##################")
error = test_trans_reordering(all_e,all_d,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   E  Fo   ##################")
error = test_trans_reordering(all_e,all_fo,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   E  Fv   ##################")
error = test_trans_reordering(all_e,all_fv,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))

print("##################   D  E    ##################")
error = test_trans_reordering(all_d,all_e,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   D  D    ##################")
error = test_trans_reordering(all_d,all_d,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   D  Fo   ##################")
error = test_trans_reordering(all_d,all_fo,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   D  Fv   ##################")
error = test_trans_reordering(all_d,all_fv,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))

print("##################   Fo E    ##################")
error = test_trans_reordering(all_fo,all_e,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   Fo D    ##################")
error = test_trans_reordering(all_fo,all_d,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   Fo Fo   ##################")
error = test_trans_reordering(all_fo,all_fo,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   Fo Fv   ##################")
error = test_trans_reordering(all_fo,all_fv,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))

print("##################   Fv E    ##################")
error = test_trans_reordering(all_fv,all_e,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   Fv D    ##################")
error = test_trans_reordering(all_fv,all_d,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   Fv Fo   ##################")
error = test_trans_reordering(all_fv,all_fo,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
print("##################   Fv Fv   ##################")
error = test_trans_reordering(all_fv,all_fv,all_occ_strings)
print("\ntotal error = {}\n\n\n\n\n\n\n\n\n\n".format(error))
