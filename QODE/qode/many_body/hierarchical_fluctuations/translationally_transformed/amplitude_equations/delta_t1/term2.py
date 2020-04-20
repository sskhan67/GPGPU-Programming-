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
            for c in range(n_virtual):
                for i in range(n_occupied):
                    for a in range(n_virtual):
                        omega_ai.amplitudes_container[n, a, i] += 1 * f_pq(a + n_occupied, n, c + n_occupied, o) * t_ai(c, o ,i ,0)

    return
