#
# header
#

# This is automatically generated code.  Do not edit directly!

from .domain_util import domain, intersection

def contraction(omega_ai, L, n_occupied, n_virtual, f_pq, v_pqrs, t_ai, t_abij):
    _domain = domain(L)
    L_0 = _domain(0)

    for n in L_0:
        L_n = _domain(n)
        for o in L_n:
            L_o = _domain(o)
            for p in intersection(L_n, L_o):
                L_p = _domain(p)
                for q in intersection(L_0, L_n, L_o, L_p):
                    L_q = _domain(q)
                    for k in range(n_occupied):
                        for c in range(n_virtual):
                            for l in range(n_occupied):
                                for i in range(n_occupied):
                                    for a in range(n_virtual):
                                        omega_ai.amplitudes_container[n, a, i] += -1 * v_pqrs(k, o, l, q, c + n_occupied, p, i, n) * t_ai(c, p ,k ,o) * t_ai(a, 0 ,l ,q)

    return
