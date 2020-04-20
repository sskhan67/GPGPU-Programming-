#import numpy as np
from .containers_classes import omega_class, fill_up_omega, single_amplitudes_class, double_amplitudes_class 



class BakerCampbellHausdorff:
    def __init__(self, limits, resources):
        self.limits    = limits
        self.resources = resources
    def computeOmega(self, H, T, textlog):
        virtuals = T.single_amplitudes.number_of_virtual_orbitals
        occupied = T.single_amplitudes.number_of_occupied_orbitals
        limit    = T.single_amplitudes.cut_off_limit
        t_ai     = single_amplitudes_class(virtuals, occupied, limit)
        t_abij   = double_amplitudes_class(virtuals, occupied, limit)
        omega = omega_class(t_ai, t_abij)
        omega = fill_up_omega(omega, H, T)
        return omega
    @staticmethod
    def Energy(omega, textlog):
        return omega.energy  # not sure about this line


