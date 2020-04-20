
#from tensors import *


#L = 1
def delta_E(L, n_virtuals, n_occupied, f_pq, v_pqrs, t_ai, t_abij):
    result = 0

    m=0


    L_0 = [i for i in range(-L, L+1)]
    L_0 = [i for i in range(-0, 0+1)]  # extra line
    domain_o = L_0

    for o in domain_o:
        for a in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o) )
            for p in domain_p:
                for i in range(n_occupied):
                    term = f_pq(i, p, a + n_occupied, o)*t_ai(a, o ,i ,p)
                    if abs(term) > 1e-6: 
                        print('product1 =', term)
                    result += 1*f_pq(i, p, a + n_occupied, o)*t_ai(a, o ,i ,p)




    L_0 = [i for i in range(-L, L+1)]
    L_0 = [i for i in range(-0, 0+1)]  # extra line
    domain_o = L_0

    for o in domain_o:
        for j in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o) )
            for p in domain_p:
                for b in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_o) )
                    for q in domain_q:
                        for a in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_o)  & set(L_q) )
                            for r in domain_r:
                                for i in range(n_occupied):
                                    term = 0.25*v_pqrs(i, r, j, o, a + n_occupied, q, b + n_occupied, p)*t_abij(a ,q  ,b ,p ,i ,r ,j ,o)
                                    if abs(term) > 1e-6: 
                                        print('product2=',term)
                                    result += 0.25*v_pqrs(i, r, j, o, a + n_occupied, q, b + n_occupied, p)*t_abij(a ,q  ,b ,p ,i ,r ,j ,o)

    print('domain_p=', domain_p)
    print('domain_q=', domain_q)
    print('domain_r=', domain_r)
    print('n_occupied =', n_occupied)
    print('n_virtual  =', n_virtuals)

    L_0 = [i for i in range(-L, L+1)]
    L_0 = [i for i in range(-0, 0+1)] # extra line
    domain_o = L_0

    for o in domain_o:
        for a in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o) )
            for p in domain_p:
                for i in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_o) )
                    for q in domain_q:
                        for j in range(n_occupied):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_o)  & set(L_q) )
                            for r in domain_r:
                                for b in range(n_virtuals):
                                    term = 0.5*v_pqrs(i, p, j, q, a + n_occupied, o, b + n_occupied, r)*t_ai(a, o ,i ,p)*t_ai(b, r ,j ,q)
                                    if abs(term) > 1e-6: 
                                        print('product3=', term)
                                    result += 0.5*v_pqrs(i, p, j, q, a + n_occupied, o, b + n_occupied, r)*t_ai(a, o ,i ,p)*t_ai(b, r ,j ,q)

    print('E_result =', result)

    return(result)
