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
from qode.fermion_field import state
from qode.many_body.coupled_cluster.attic.operator_commutation import HBar_on_HF



class NestedCommutators(object):
	def __init__(self, cisd_strings, resources):
		self.cisd_strings = cisd_strings
		self.resources    = resources
	def __call__(self, H, T, textlog):
		omega0_state = HBar_on_HF(H, T, self.cisd_strings, self.resources, textlog)
		HBar = deepcopy(T)
		HBar.t_amp_obj = omega0_state
		return HBar



class PseudoEnergy(object):
	def __init__(self, cisd_strings, resources):
		self.cisd_strings = cisd_strings
		self.resources    = resources
	def __call__(self, HBar, textlog):
		return state.dot(state.state(self.cisd_strings,0), HBar.t_amp_obj)



class BakerCampbellHausdorff(object):
	def __init__(self, cisd_strings, resources):
		self.BCH = NestedCommutators(cisd_strings, resources)
		self.E   = PseudoEnergy(cisd_strings, resources)



