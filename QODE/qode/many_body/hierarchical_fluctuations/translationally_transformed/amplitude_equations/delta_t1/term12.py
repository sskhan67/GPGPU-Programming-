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
        for o in L_0:
            L_o = _domain(o)
            for p in L_o:
                L_p = _domain(p)
                for q in intersection(L_0, L_n, L_o, L_p):
                    L_q = _domain(q)
                    for r in intersection(L_0, L_n, L_o, L_p, L_q):
                        L_r = _domain(r)
                        for k in range(n_occupied):
                            for c in range(n_virtual):
                                for d in range(n_virtual):
                                    for l in range(n_occupied):
                                        for i in range(n_occupied):
                                            for a in range(n_virtual):
                                                omega_ai.amplitudes_container[n, a, i] += 1 * v_pqrs(k, o, l, r, c + n_occupied, p, d + n_occupied, q) * t_ai(c, p ,k ,o) * t_abij(d ,q  ,a ,n ,l ,r ,i ,0)

    return
