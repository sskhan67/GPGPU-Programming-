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

the first step is to run ../main.py with the dump flag un-disabled in qode/many_body/nested_operator/baker_campbell_hausdorff/bch.py (and using direct commutation).
also dump the delT's in ../hamiltonian.py.

copy the energy from the 4th cycle (the dumped one) into reference/energy.txt and move the pickled files into the reference/ directory (copy the 4th delT as reference/delT_raw.pickled).
double-check with diff_ops.py that H.picked and X_0.pickled hold equivalent operators.

then with dump and direct commutation still on, run build_Omega.py ... it should print the same energy three times and the error against the reference Omega and delT operators.
then can turn off dumps (X's and delT's) and direct commutation.
can also manually diff_ops.py the newly dumped files against the reference files. ... then I just deleted the newly dumped files.
these check the integrity of H, T and the procedure to build omega.

in order to check the integrity of the procedure that will be used to build the decomposed debugging data, run check_decomp.py on each of the pickled operator files now in reference/.
It should check that decomposing and resumming is the same as the original and print the (hopefully small) error.

copy H.pickled T.pickled and Omega.pickled from reference/ to masked/
copy each of the reference/X_<n>.pickled files to masked/X_<n>/total.pickled

run decompose_commutator.py with the arguments 1, 2, 3, and 4 to generate from*.pickled files in each of the masked/X_<n>/ directories.
These files give the parts of the commutators that are made by commuting the respective blocks of X_<n-1> with T.
It also generates a check.pickled file in each directory, which is not decomposed, but uses othewise the same means, so it should
be identical to the respective total.pickled, which can be checked using diff_ops.py.

Finally, sum.py, called with the arguments 1, 2, 3 and 4 will sum up all of the from*.pickled files in a given directory
and make the file sum.pickled there, which should represent the same operator as in total.pickled and check.pickled in the respective
directory.  The errors between sum.picked and total.pickled will be larger than between check.pickled and total.pickled
because check and total were ultimately generated doing all the flops in the same order (exactly, I think, but only with
some extra dumps to disk in between).  In my experience, they were substantially larger, but still around 1e-16 or less.  I think
this is basically due to the roundoff precision of the machine.  Even if some smaller operators were generated, there
may have been cancelling terms in the process, and so the error is approximately the roundoff threshold.
Even in the worst case, the error was still 10 orders of magnitude smaller than the norm.

If all the checking above passes, then one can believe that the from*.pickled files are sound and use those for debugging other implementations
of operator commutation.

