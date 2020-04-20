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
        for o in intersection(L_0, L_n):
            L_o = _domain(o)
            for p in intersection(L_0, L_n, L_o):
                L_p = _domain(p)
                for q in intersection(L_0, L_n, L_o, L_p):
                    L_q = _domain(q)
                    for d in range(n_virtual):
                        for k in range(n_occupied):
                            for c in range(n_virtual):
                                for i in range(n_occupied):
                                    for a in range(n_virtual):
                                        omega_ai.amplitudes_container[n, a, i] += 0.5 * v_pqrs(k, p, a + n_occupied, n, c + n_occupied, q, d + n_occupied, o) * t_abij(c ,q  ,d ,o ,k ,p ,i ,0)

    return
