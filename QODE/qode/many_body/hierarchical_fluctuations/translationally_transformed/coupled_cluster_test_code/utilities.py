import numpy
from containers_classes import single_amplitudes_class, double_amplitudes_class, amplitudes_class
from containers_classes import one_body_hamiltonian_class, two_body_hamiltonian_class, hamiltonian_class
from containers_classes import omega_class, fill_up_omega, delta_t, fill_up_one_body_hamiltonian
from containers_classes import fill_up_two_body_hamiltonian,fill_up_single_amplitudes, get_omega
from containers_classes import fill_up_double_amplitudes,inv_F_singles, inv_F_doubles, inv_F_class, fill_up_inv_F
from operations_classes import dot, add_to, scale, copy, act_on_vec


class space_traits_class(object):
    def __init__(self):
        self.field = numpy.float64
    @staticmethod
    def dot(v,w):
        return dot(v,w)
    @staticmethod
    def add_to(v, w, n=1.0):
        add_to(v, w, n)
    @staticmethod
    def scale(n, v):
        scale(n, v)
    @staticmethod
    def copy(v):
        u = copy(v)  
        return u
    @staticmethod
    def act_on_vec(op, v):
        return act_on_vec(op, v)   # This is a single tasker, only here for acting orb energy denominators
        ############################# not needed for this #################################
    def check_member(self,v):        pass
    def check_lin_op(self,op):       return False
    @staticmethod
    def function_on_diags(func,op):  raise NotImplementedError
    @staticmethod
    def back_act_on_vec(v,op):       raise NotImplementedError
    @staticmethod
    def diagonal(op):                raise NotImplementedError
space_traits = space_traits_class()

