import numpy as np
from qode.many_body.coupled_cluster.ccsd import CCSD  
from qode.util import parallel, output, textlog, indent

from qode.many_body.hierarchical_fluctuations.translationally_transformed.baker_campbell_hausdorff import BakerCampbellHausdorff
from qode.many_body.hierarchical_fluctuations.translationally_transformed.utilities import space_traits
from qode.many_body.hierarchical_fluctuations.translationally_transformed.containers_classes import single_amplitudes_class, double_amplitudes_class, amplitudes_class
from qode.many_body.hierarchical_fluctuations.translationally_transformed.containers_classes import one_body_hamiltonian_class, two_body_hamiltonian_class, hamiltonian_class
from qode.many_body.hierarchical_fluctuations.translationally_transformed.containers_classes import omega_class, fill_up_omega, fill_up_one_body_hamiltonian
from qode.many_body.hierarchical_fluctuations.translationally_transformed.containers_classes import fill_up_two_body_hamiltonian, fill_up_hamiltonian
from qode.many_body.hierarchical_fluctuations.translationally_transformed.containers_classes import inv_F_singles, inv_F_doubles, inv_F_class, fill_up_inv_F
from qode.many_body.hierarchical_fluctuations.translationally_transformed.containers_classes import load_parameters

# setting up parameters
distance = 3.0
cut_off_limit = 0
electrons = 2
parameters = load_parameters('helium', 'orthogonal', '6-31G', cut_off_limit, distance, 'semiMO')

# deducing extra paramters from inputted values
number_of_occupied_orbitals = electrons
number_of_virtual_orbitals  = int(parameters[0][0][0].shape[0]*2 - number_of_occupied_orbitals) 

# creating hamiltonian class
f_pq = one_body_hamiltonian_class(number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit)
v_pqrs = two_body_hamiltonian_class(number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit)
H = hamiltonian_class(f_pq, v_pqrs)

# creating amplitudes class
t_ai = single_amplitudes_class(number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit)
t_abij = double_amplitudes_class(number_of_virtual_orbitals, number_of_occupied_orbitals, cut_off_limit)
T = amplitudes_class(t_ai, t_abij)

# creating denominators class
f1 = inv_F_singles(number_of_occupied_orbitals, number_of_virtual_orbitals, cut_off_limit)
f2 = inv_F_doubles(number_of_occupied_orbitals, number_of_virtual_orbitals, cut_off_limit)
inv_F = inv_F_class(f1, f2)

# filling up classes
inv_F = fill_up_inv_F(inv_F, parameters)
H = fill_up_hamiltonian(H, parameters)

# amplitude optimization
out, resources = output(log=textlog(echo=True)), parallel.resources(1)   # this line was taken from Yuhong's code
BCH = BakerCampbellHausdorff(cut_off_limit, resources)
E = CCSD(H, T, inv_F, BCH, space_traits, out.log.sublog(), thresh=1e-6)
print( E )
