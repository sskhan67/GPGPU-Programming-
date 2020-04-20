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
In the original/ folder:

This is the code that Yuhong used to produce the timing data for recoupled oscillator molecules for the first paper,
but, which unfortunately used the buggy local_operator (autoQode) code.  The only thing which has been changed is
the paths to the modules he used, since those moved.  There is the subtle detail that the local_operator code
had one bug fixed that lies outside of autoQode.  In the space_traits class for the cluster operator, the add_to
function was scaling the incremented vector in place before adding it to the target.  However, it has been verified
that this made no difference in this application (I had some reference data from before the bug fix; every energy
agreed to every digit).  I guess this is because increments were discarded after their only use
(makes sense for omega, but what about the insides of DIIS)? ... regardless, a dead bug is better than a live one, 
and the fact that this changed does not affect the logic below.

In the validate/ folder there are two subfolders:

	old/
The code in here is nearly identical to that in the original/ folder above.  The only difference is that the outer-most
loops have been modified to only run a substantial but small subset of the total number of calculations.
	new/
The code in this directory is nearly identical to that in the top Applications/recoupled-osc-ccsd/ folder except that
it uses the (buggy) local_operator code and it only runs the same subset of the total number of timing jobs.
Although this code is much more cleaned up and restructured (making it easier to swap out the local_operator code
for its replacement), in principle, it does the same thing as the code in the old/ folder.
By directly diff-ing outputs of code from old/ and new/, it is seen that timings are similar and that any energies
which differ do so at the machine threshold.  Therefore, the structure of the new code has been validated.

In the reproduce/ folder:

This code is nearly identical to that in validate/new, except that it runs all calculations, which makes it nearly
identical to the that in the top Applications/recoupled-osc-ccsd/ folder except that it uses the (buggy) local_operator
code.  The purpose of this code was to, just for good measure, completely duplicate Yuhong's work, including the
final figure that results, just to make sure I am doing absolutely everything right, including parsing of the data,
before swapping local_operator for its replacement.

The result was affirmative, aside from the fact that (hardware-related?) crashes caused me to give up after 22 molecules
for the primitive variant.  For the excitonic plot this validation result is visually the same as the prior plot, but,
most importantly on the visual front, both the excitonic and primitive results fall relatively nicely on single-monomial
best-fit curves.  The results from this validation run is that these curves are 0.035sec*N**2.31 and 0.046sec*N**3.00 for the
excitonic and primitive curves, respectively, whereas the prior results were 0.036sec*N**2.30 and 0.054sec*N**2.92,
which is in very good agreement (and meaningful, since the plots look good), especially considering that the fits use data
from only systems with 15 or more molecules.  So the validation result should be slightly different for the primitive case,
both due to substantially worse statistics and the fact that I know there are overhead (0th-order) terms that make the precise 
result dependent on which data range is selected (which is why I previously rejected smaller calcs).
