#
# header
#

# This is automatically generated code.  Do not edit directly!

import ctypes
from qode.util import get_ctype
from .location import QodeHome	# DIRTY and DANGEROUS: fix this!
from .domain_util import domain, intersection

def contraction(omega_ai, L, n_occupied, n_virtual, f_pq, v_pqrs, t_ai, t_abij):
    # gcc -std=c99 -c -fPIC -O2 term14_inner.c -o term14_inner.o
    # gcc -shared -O2 term14_inner.o -o term14_inner.so
    term14_inner = ctypes.cdll.LoadLibrary(QodeHome + "/qode/many_body/hierarchical_fluctuations/translationally_transformed/amplitude_equations/delta_t1/term14_inner.so")
    C_n_occupied = get_ctype(n_occupied, "int")
    C_n_virtual  = get_ctype(n_virtual,  "int")
    C_omega_ai   = get_ctype(omega_ai.raw_storage(), "double*")
    C_f_pq       = get_ctype(    f_pq.raw_storage(), "double*")
    C_v_pqrs     = get_ctype(  v_pqrs.raw_storage(), "double*")
    C_t_ai       = get_ctype(    t_ai.raw_storage(), "double*")
    C_t_abij     = get_ctype(  t_abij.raw_storage(), "double*")

    _domain = domain(L)
    L_0 = _domain(0)

    for n in L_0:
        L_n = _domain(n)
        for o in L_n:
            L_o = _domain(o)
            for p in intersection(L_n, L_o):
                L_p = _domain(p)
                for q in intersection(L_n, L_o, L_p):
                    L_q = _domain(q)
                    for r in intersection(L_0, L_o, L_p, L_q):
                        L_r = _domain(r)
                        #for k in range(n_occupied):
                        #    for c in range(n_virtual):
                        #        for l in range(n_occupied):
                        #            for d in range(n_virtual):
                        #                for i in range(n_occupied):
                        #                    for a in range(n_virtual):
                        #                        omega_ai.amplitudes_container[n, a, i] += -0.5 * v_pqrs(k, o, l, q, c + n_occupied, p, d + n_occupied, r) * t_abij(c ,p  ,a ,n ,k ,o ,l ,q) * t_ai(d, r ,i ,0)
                        omega_ai_offset = omega_ai.raw_offset(n, 0)         # must be saved to a local variable or else not there when C goes to look for it
                        f_pq_offset     = -1                                # must be saved to a local variable or else not there when C goes to look for it
                        v_pqrs_offset   =   v_pqrs.raw_offset(o, q, p, r)   # must be saved to a local variable or else not there when C goes to look for it
                        t_ai_offset     =     t_ai.raw_offset(r, 0)         # must be saved to a local variable or else not there when C goes to look for it
                        t_abij_offset   =   t_abij.raw_offset(p, n, o, q)   # must be saved to a local variable or else not there when C goes to look for it
                        term14_inner.loops(C_n_occupied, C_n_virtual, C_omega_ai, C_f_pq, C_v_pqrs, C_t_ai, C_t_abij, get_ctype(omega_ai_offset, "int"), get_ctype(f_pq_offset, "int"), get_ctype(v_pqrs_offset, "int"), get_ctype(t_ai_offset, "int"), get_ctype(t_abij_offset, "int"))

    return
