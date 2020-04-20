

def delta_t2(i,j,a,b,n_1,n_2,n_3,L, n_virtuals, n_occupied, f_pq, v_pqrs, t_ai, t_abij):
    result = 0

    print(t_abij.amplitudes_container[0,0], 'slice of t2')
    m=0
    previous_result = 0

    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]


    result += 1*v_pqrs(a + n_occupied, n_2, b + n_occupied, n_3, i, m, j, n_1) # 1
    term_1 = result - previous_result
    previous_result = result
    print('term1=', term_1)


    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            result += 1*f_pq(b + n_occupied, n_3, c + n_occupied, o)*t_abij(a, n_2 , c, o, i, m, j, n_1)  # 2
    term_2 = result - previous_result
    previous_result = result
    print('term2=', term_2)


    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            result += -1*f_pq(a + n_occupied, n_2, c + n_occupied, o)*t_abij(b, n_3 , c, o, i, m, j, n_1)  # 3
    term_3 = result - previous_result
    previous_result = result
    print('term3=', term_3)


    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            result += -1*f_pq(k, o, j, n_1)*t_abij(a, n_2 , b, n_3, i, m, k, o)  # 4
    term_4 = result - previous_result
    previous_result = result
    print('term4=', term_4)


    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            result += 1*f_pq(k, o, i, m)*t_abij(a, n_2 , b, n_3, j, n_1, k, o)  # 5
    term_5 = result - previous_result
    previous_result = result
    print('term5=', term_5)


    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for l in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    result += 0.5*v_pqrs(k, p, l, o, i, m, j, n_1)*t_abij(a, n_2 , b, n_3, k, p, l, o)  # 6
    term_6 = result - previous_result
    previous_result = result
    print('term6=', term_6)


    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_n_2) )
            for p in domain_p:
                for d in range(n_virtuals):
                    result += 0.5*v_pqrs(a + n_occupied, n_2, b + n_occupied, n_3, c + n_occupied, o, d + n_occupied, p)*t_abij(c, o , d, p, i, m, j, n_1)  # 7
    term_7 = result - previous_result
    previous_result = result
    print('term7=', term_7)


    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    print('-----------------------------------------')
    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    if v_pqrs(k, p, b + n_occupied, n_3, c + n_occupied, o, j, n_1) != 0.0 and t_abij(a, n_2 , c, o, i, m, k, p) != 0.0:
                        print('v=', v_pqrs(k, p, b + n_occupied, n_3, c + n_occupied, o, j, n_1), 't=', t_abij(a, n_2 , c, o, i, m, k, p), 'multiplied in loop')
                        print('a=', a, 'c=', c, 'i=', i, 'k=', k, 'n_2=', n_2, 'o=', o, 'm=', m, 'p=', p)
                    result += 1*v_pqrs(k, p, b + n_occupied, n_3, c + n_occupied, o, j, n_1)*t_abij(a, n_2 , c, o, i, m, k, p)  # 8
    term_8 = result - previous_result
    previous_result = result
    print('term8=', term_8)
    print('-----------------------------------------')

    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    result += -1*v_pqrs(k, p, b + n_occupied, n_3, c + n_occupied, o, i, m)*t_abij(a, n_2 , c, o, j, n_1, k, p)  # 9
    term_9 = result - previous_result
    previous_result = result
    print('term9=', term_9)


    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    result += -1*v_pqrs(k, p, a + n_occupied, n_2, c + n_occupied, o, j, n_1)*t_abij(b, n_3 , c, o, i, m, k, p)  # 10
    term_10 = result - previous_result
    previous_result = result
    print('term10=', term_10)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    result += 1*v_pqrs(k, p, a + n_occupied, n_2, c + n_occupied, o, i, m)*t_abij(b, n_3 , c, o, j, n_1, k, p)  #11
    term_11 = result - previous_result
    previous_result = result
    print('term11', term_11)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            result += 1*v_pqrs(a + n_occupied, n_2, b + n_occupied, n_3, c + n_occupied, o, j, n_1)*t_ai(c, o, i, m)  #12
    term_12 = result - previous_result
    previous_result = result
    print('term12', term_12)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            result += -1*v_pqrs(a + n_occupied, n_2, b + n_occupied, n_3, c + n_occupied, o, i, m)*t_ai(c, o, j, n_1)  #13
    term_13 = result - previous_result
    previous_result = result
    print('term13', term_13)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            result += -1*v_pqrs(k, o, b + n_occupied, n_3, i, m, j, n_1)*t_ai(a, n_2, k, o)  #14
    term_14 = result - previous_result
    previous_result = result
    print('term14', term_14)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            result += 1*v_pqrs(k, o, a + n_occupied, n_2, i, m, j, n_1)*t_ai(b, n_3, k, o)  #15
    term_15 = result - previous_result
    previous_result = result
    print('term15', term_15)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_q)  & set(L_p) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += 0.5*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_abij(a, n_2 , c, o, i, m, k, p)*t_abij(d, q , b, n_3, l, r, j, n_1)  #16
    term_16 = result - previous_result
    previous_result = result
    print('term16', term_16)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_3) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_1)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_n_2) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += -0.5*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_abij(b, n_3 , c, o, i, m, k, p)*t_abij(d, q , a, n_2, l, r, j, n_1)  #17
    term_17 = result - previous_result
    previous_result = result
    print('term17', term_17)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_3)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_3)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_m) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += -0.5*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_abij(a, n_2 , c, o, j, n_1, k, p)*t_abij(d, q , b, n_3, l, r, i, m)  #18
    term_18 = result - previous_result
    previous_result = result
    print('term18', term_18)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_3) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_m)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_2)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_m) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += 0.5*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_abij(b, n_3 , c, o, j, n_1, k, p)*t_abij(d, q , a, n_2, l, r, i, m)  #19
    term_19 = result - previous_result
    previous_result = result
    print('term19', term_19)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_o) )
            for p in domain_p:
                for d in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_3)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_n_2) )
                            for r in domain_r:
                                for k in range(n_occupied):
                                    result += 0.25*v_pqrs(k, r, l, q, c + n_occupied, o, d + n_occupied, p)*t_abij(c, o , d, p, i, m, j, n_1)*t_abij(a, n_2 , b, n_3, k, r, l, q)  #20
    term_20 = result - previous_result
    previous_result = result
    print('term20', term_20)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for d in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_3)  & set(L_o) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_n_3)  & set(L_q)  & set(L_o) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += -0.5*v_pqrs(k, q, l, r, c + n_occupied, o, d + n_occupied, p)*t_abij(a, n_2 , c, o, i, m, j, n_1)*t_abij(b, n_3 , d, p, k, q, l, r)  #21
    term_21 = result - previous_result
    previous_result = result
    print('term21', term_21)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for d in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_q)  & set(L_n_2)  & set(L_o) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += 0.5*v_pqrs(k, q, l, r, c + n_occupied, o, d + n_occupied, p)*t_abij(b, n_3 , c, o, i, m, j, n_1)*t_abij(a, n_2 , d, p, k, q, l, r)  #22
    term_22 = result - previous_result
    previous_result = result
    print('term22', term_22)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_m)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_o) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_n_1)  & set(L_q)  & set(L_o) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += -0.5*v_pqrs(k, o, l, r, c + n_occupied, p, d + n_occupied, q)*t_abij(a, n_2 , b, n_3, i, m, k, o)*t_abij(c, p , d, q, j, n_1, l, r)  #23
    term_23 = result - previous_result
    previous_result = result
    print('term23', term_23)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_o) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_m)  & set(L_q)  & set(L_o) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += 0.5*v_pqrs(k, o, l, r, c + n_occupied, p, d + n_occupied, q)*t_abij(a, n_2 , b, n_3, j, n_1, k, o)*t_abij(c, p , d, q, i, m, l, r)  #24
    term_24 = result - previous_result
    previous_result = result
    print('term24', term_24)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for l in range(n_occupied):
                    result += 0.5*v_pqrs(k, o, l, p, i, m, j, n_1)*t_ai(a, n_2, k, o)*t_ai(b, n_3, l, p)  #25
    term_25 = result - previous_result
    previous_result = result
    print('term25', term_25)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1)  & set(L_n_3) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for l in range(n_occupied):
                    result += -0.5*v_pqrs(k, o, l, p, i, m, j, n_1)*t_ai(b, n_3, k, o)*t_ai(a, n_2, l, p)  #26
    term_26 = result - previous_result
    previous_result = result
    print('term26', term_26)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_m)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for d in range(n_virtuals):
                    result += 0.5*v_pqrs(a + n_occupied, n_2, b + n_occupied, n_3, c + n_occupied, o, d + n_occupied, p)*t_ai(c, o, i, m)*t_ai(d, p, j, n_1)  #27
    term_27 = result - previous_result
    previous_result = result
    print('term27', term_27)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for d in range(n_virtuals):
                    result += -0.5*v_pqrs(a + n_occupied, n_2, b + n_occupied, n_3, c + n_occupied, o, d + n_occupied, p)*t_ai(c, o, j, n_1)*t_ai(d, p, i, m)  #28
    term_28 = result - previous_result
    previous_result = result
    print('term28', term_28)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_m)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for c in range(n_virtuals):
                    result += -1*v_pqrs(k, o, b + n_occupied, n_3, i, m, c + n_occupied, p)*t_ai(a, n_2, k, o)*t_ai(c, p, j, n_1)  #29
    term_29 = result - previous_result
    previous_result = result
    print('term29', term_29)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_m)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    result += 1*v_pqrs(k, o, a + n_occupied, n_2, i, m, c + n_occupied, p)*t_ai(b, n_3, k, o)*t_ai(c, p, j, n_1)  #30
    term_30 = result - previous_result
    previous_result = result
    print('term30', term_30)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for c in range(n_virtuals):
                    result += 1*v_pqrs(k, o, b + n_occupied, n_3, j, n_1, c + n_occupied, p)*t_ai(a, n_2, k, o)*t_ai(c, p, i, m)  #31
    term_31 = result - previous_result
    previous_result = result
    print('term31', term_31)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    result += -1*v_pqrs(k, o, a + n_occupied, n_2, j, n_1, c + n_occupied, p)*t_ai(b, n_3, k, o)*t_ai(c, p, i, m)  #32
    term_32 = result - previous_result
    previous_result = result
    print('term32', term_32)




    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = L_n_2

    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for c in range(n_virtuals):
                    result += 1*f_pq(k, o, c + n_occupied, p)*t_ai(a, n_2, k, o)*t_abij(b, n_3 , c, p, i, m, j, n_1)  #33
    term_33 = result - previous_result
    previous_result = result
    print('term33', term_33)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_n_3

    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    result += -1*f_pq(k, o, c + n_occupied, p)*t_ai(b, n_3, k, o)*t_abij(a, n_2 , c, p, i, m, j, n_1)  #34
    term_34 = result - previous_result
    previous_result = result
    print('term34', term_34)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_m

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    result += 1*f_pq(k, p, c + n_occupied, o)*t_ai(c, o, i, m)*t_abij(a, n_2 , b, n_3, j, n_1, k, p)  #35
    term_35 = result - previous_result
    previous_result = result
    print('term35', term_35)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_n_1

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    result += -1*f_pq(k, p, c + n_occupied, o)*t_ai(c, o, j, n_1)*t_abij(a, n_2 , b, n_3, i, m, k, p)  #36
    term_36 = result - previous_result
    previous_result = result
    print('term36', term_36)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_m

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_p)  & set(L_n_2) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            result += -1*v_pqrs(k, p, l, q, c + n_occupied, o, i, m)*t_ai(c, o, k, p)*t_abij(a, n_2 , b, n_3, l, q, j, n_1)  #37
    term_37 = result - previous_result
    previous_result = result
    print('term37', term_37)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_n_1

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_p)  & set(L_n_2) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            result += 1*v_pqrs(k, p, l, q, c + n_occupied, o, j, n_1)*t_ai(c, o, k, p)*t_abij(a, n_2 , b, n_3, l, q, i, m)  #38
    term_38 = result - previous_result
    previous_result = result
    print('term38', term_38)




    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = L_n_2

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_p)  & set(L_n_2) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            result += 1*v_pqrs(k, p, a + n_occupied, n_2, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, k, p)*t_abij(d, q , b, n_3, i, m, j, n_1)  #39
    term_39 = result - previous_result
    previous_result = result
    print('term39', term_39)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_n_3

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_p)  & set(L_n_2) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            result += -1*v_pqrs(k, p, b + n_occupied, n_3, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, k, p)*t_abij(d, q , a, n_2, i, m, j, n_1)  #40
    term_40 = result - previous_result
    previous_result = result
    print('term40', term_40)




    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_2) )


    for o in domain_o:
        for d in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_p)  & set(L_n_2) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += 1*v_pqrs(a + n_occupied, n_2, k, q, d + n_occupied, o, c + n_occupied, p)*t_ai(d, o, i, m)*t_abij(b, n_3 , c, p, j, n_1, k, q)  #41
    term_41 = result - previous_result
    previous_result = result
    print('term41', term_41)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_3) )


    for o in domain_o:
        for d in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_p)  & set(L_n_2) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += -1*v_pqrs(b + n_occupied, n_3, k, q, d + n_occupied, o, c + n_occupied, p)*t_ai(d, o, i, m)*t_abij(a, n_2 , c, p, j, n_1, k, q)  #42
    term_42 = result - previous_result
    previous_result = result
    print('term42', term_42)




    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_2) )


    for o in domain_o:
        for d in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_3)  & set(L_n_2)  & set(L_o)  & set(L_p)  & set(L_m) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += -1*v_pqrs(a + n_occupied, n_2, k, q, d + n_occupied, o, c + n_occupied, p)*t_ai(d, o, j, n_1)*t_abij(b, n_3 , c, p, i, m, k, q)  #43
    term_43 = result - previous_result
    previous_result = result
    print('term43', term_43)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_3) )


    for o in domain_o:
        for d in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_3)  & set(L_n_2)  & set(L_o)  & set(L_p)  & set(L_m) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += 1*v_pqrs(b + n_occupied, n_3, k, q, d + n_occupied, o, c + n_occupied, p)*t_ai(d, o, j, n_1)*t_abij(a, n_2 , c, p, i, m, k, q)  #44
    term_44 = result - previous_result
    previous_result = result
    print('term44', term_44)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_2) )


    for o in domain_o:
        for l in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_p)  & set(L_m) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += 1*v_pqrs(k, q, l, o, i, m, c + n_occupied, p)*t_ai(a, n_2, l, o)*t_abij(b, n_3 , c, p, j, n_1, k, q)  #45
    term_45 = result - previous_result
    previous_result = result
    print('term45', term_45)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_3) )


    for o in domain_o:
        for l in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_2)  & set(L_o)  & set(L_p)  & set(L_m) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += -1*v_pqrs(k, q, l, o, i, m, c + n_occupied, p)*t_ai(b, n_3, l, o)*t_abij(a, n_2 , c, p, j, n_1, k, q)  #46
    term_46 = result - previous_result
    previous_result = result
    print('term46', term_46)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_2) )


    for o in domain_o:
        for l in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_p)  & set(L_m) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += -1*v_pqrs(k, q, l, o, j, n_1, c + n_occupied, p)*t_ai(a, n_2, l, o)*t_abij(b, n_3 , c, p, i, m, k, q)  #47
    term_47 = result - previous_result
    previous_result = result
    print('term47', term_47)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_3) )


    for o in domain_o:
        for l in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_2)  & set(L_o)  & set(L_p)  & set(L_m) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += 1*v_pqrs(k, q, l, o, j, n_1, c + n_occupied, p)*t_ai(b, n_3, l, o)*t_abij(a, n_2 , c, p, i, m, k, q)  #48
    term_48 = result - previous_result
    previous_result = result
    print('term48', term_48)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for l in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_p)  & set(L_n_2) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += 0.5*v_pqrs(k, q, l, p, c + n_occupied, o, j, n_1)*t_ai(c, o, i, m)*t_abij(a, n_2 , b, n_3, k, q, l, p)  #49
    term_49 = result - previous_result
    previous_result = result
    print('term49', term_49)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for l in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_3)  & set(L_n_2)  & set(L_o)  & set(L_p)  & set(L_m) )
                    for q in domain_q:
                        for k in range(n_occupied):
                            result += -0.5*v_pqrs(k, q, l, p, c + n_occupied, o, i, m)*t_ai(c, o, j, n_1)*t_abij(a, n_2 , b, n_3, k, q, l, p)  #50
    term_50 = result - previous_result
    previous_result = result
    print('term50', term_50)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_p)  & set(L_m) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            result += -0.5*v_pqrs(k, o, b + n_occupied, n_3, c + n_occupied, p, d + n_occupied, q)*t_ai(a, n_2, k, o)*t_abij(c, p , d, q, i, m, j, n_1)  #51
    term_51 = result - previous_result
    previous_result = result
    print('term51', term_51)




    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_n_3)  & set(L_n_2) )


    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for c in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_n_1)  & set(L_n_2)  & set(L_o)  & set(L_p)  & set(L_m) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            result += 0.5*v_pqrs(k, o, a + n_occupied, n_2, c + n_occupied, p, d + n_occupied, q)*t_ai(b, n_3, k, o)*t_abij(c, p , d, q, i, m, j, n_1)  #52
    term_52 = result - previous_result
    previous_result = result
    print('term52', term_52)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_m)  & set(L_n_3) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            result += -0.5*v_pqrs(k, p, b + n_occupied, n_3, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, i, m)*t_ai(a, n_2, k, p)*t_ai(d, q, j, n_1)  #53
    term_53 = result - previous_result
    previous_result = result
    print('term53', term_53)




    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = list( set(L_m)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            result += 0.5*v_pqrs(k, p, a + n_occupied, n_2, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, i, m)*t_ai(b, n_3, k, p)*t_ai(d, q, j, n_1)  #54
    term_54 = result - previous_result
    previous_result = result
    print('term54', term_54)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_3) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_3)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            result += 0.5*v_pqrs(k, p, b + n_occupied, n_3, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, j, n_1)*t_ai(a, n_2, k, p)*t_ai(d, q, i, m)  #55
    term_55 = result - previous_result
    previous_result = result
    print('term55', term_55)




    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    domain_o = list( set(L_n_1)  & set(L_n_2) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_m)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            result += -0.5*v_pqrs(k, p, a + n_occupied, n_2, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, j, n_1)*t_ai(b, n_3, k, p)*t_ai(d, q, i, m)  #56
    term_56 = result - previous_result
    previous_result = result
    print('term56', term_56)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            result += 0.5*v_pqrs(k, p, l, q, c + n_occupied, o, j, n_1)*t_ai(c, o, i, m)*t_ai(a, n_2, k, p)*t_ai(b, n_3, l, q)  #57
    term_57 = result - previous_result
    previous_result = result
    print('term57', term_57)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            result += -0.5*v_pqrs(k, p, l, q, c + n_occupied, o, j, n_1)*t_ai(c, o, i, m)*t_ai(b, n_3, k, p)*t_ai(a, n_2, l, q)  #58
    term_58 = result - previous_result
    previous_result = result
    print('term58', term_58)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_3)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            result += -0.5*v_pqrs(k, p, l, q, c + n_occupied, o, i, m)*t_ai(c, o, j, n_1)*t_ai(a, n_2, k, p)*t_ai(b, n_3, l, q)  #59
    term_59 = result - previous_result
    previous_result = result
    print('term59', term_59)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = list( set(L_m)  & set(L_n_1) )


    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_m)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            result += 0.5*v_pqrs(k, p, l, q, c + n_occupied, o, i, m)*t_ai(c, o, j, n_1)*t_ai(b, n_3, k, p)*t_ai(a, n_2, l, q)  #60
    term_60 = result - previous_result
    previous_result = result
    print('term60', term_60)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_0 = [i for i in range(-L, L+1)]
    domain_o = L_0

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_n_2) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += -1*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, k, p)*t_ai(d, q, i, m)*t_abij(a, n_2 , b, n_3, l, r, j, n_1)  #61
    term_61 = result - previous_result
    previous_result = result
    print('term61', term_61)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_0 = [i for i in range(-L, L+1)]
    domain_o = L_0

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_3)  & set(L_o)  & set(L_m)  & set(L_q)  & set(L_p)  & set(L_n_2) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += 1*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, k, p)*t_ai(d, q, j, n_1)*t_abij(a, n_2 , b, n_3, l, r, i, m)  #62
    term_62 = result - previous_result
    previous_result = result
    print('term62', term_62)




    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_0 = [i for i in range(-L, L+1)]
    domain_o = L_0

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_m) )
                            for r in domain_r:
                                for d in range(n_virtuals):
                                    result += -1*v_pqrs(k, p, l, q, c + n_occupied, o, d + n_occupied, r)*t_ai(c, o, k, p)*t_ai(a, n_2, l, q)*t_abij(d, r , b, n_3, i, m, j, n_1)  #63
    term_63 = result - previous_result
    previous_result = result
    print('term63', term_63)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_0 = [i for i in range(-L, L+1)]
    domain_o = L_0

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_3)  & set(L_o) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_1)  & set(L_o)  & set(L_m)  & set(L_q)  & set(L_p)  & set(L_n_2) )
                            for r in domain_r:
                                for d in range(n_virtuals):
                                    result += 1*v_pqrs(k, p, l, q, c + n_occupied, o, d + n_occupied, r)*t_ai(c, o, k, p)*t_ai(b, n_3, l, q)*t_abij(d, r , a, n_2, i, m, j, n_1)  #64
    term_64 = result - previous_result
    previous_result = result
    print('term64', term_64)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_m

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_1)  & set(L_o) )
            for p in domain_p:
                for d in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_3)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_n_2) )
                            for r in domain_r:
                                for k in range(n_occupied):
                                    result += 0.25*v_pqrs(k, r, l, q, c + n_occupied, o, d + n_occupied, p)*t_ai(c, o, i, m)*t_ai(d, p, j, n_1)*t_abij(a, n_2 , b, n_3, k, r, l, q)  #65
    term_65 = result - previous_result
    previous_result = result
    print('term65', term_65)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_n_1

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_m)  & set(L_o) )
            for p in domain_p:
                for d in range(n_virtuals):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_3)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for l in range(n_occupied):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_3)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_n_2) )
                            for r in domain_r:
                                for k in range(n_occupied):
                                    result += -0.25*v_pqrs(k, r, l, q, c + n_occupied, o, d + n_occupied, p)*t_ai(c, o, j, n_1)*t_ai(d, p, i, m)*t_abij(a, n_2 , b, n_3, k, r, l, q)  #66
    term_66 = result - previous_result
    previous_result = result
    print('term66', term_66)




    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = L_n_2

    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for l in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for c in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_1)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_m) )
                            for r in domain_r:
                                for d in range(n_virtuals):
                                    result += 0.25*v_pqrs(k, o, l, p, c + n_occupied, q, d + n_occupied, r)*t_ai(a, n_2, k, o)*t_ai(b, n_3, l, p)*t_abij(c, q , d, r, i, m, j, n_1)  #67
    term_67 = result - previous_result
    previous_result = result
    print('term67', term_67)




    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    domain_o = L_n_3

    for o in domain_o:
        for k in range(n_occupied):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for l in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for c in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_1)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_m) )
                            for r in domain_r:
                                for d in range(n_virtuals):
                                    result += -0.25*v_pqrs(k, o, l, p, c + n_occupied, q, d + n_occupied, r)*t_ai(b, n_3, k, o)*t_ai(a, n_2, l, p)*t_abij(c, q , d, r, i, m, j, n_1)  #68
    term_68 = result - previous_result
    previous_result = result
    print('term68', term_68)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_m

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for l in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_1)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_n_2) )
                            for r in domain_r:
                                for k in range(n_occupied):
                                    result += 1*v_pqrs(k, r, l, p, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, i, m)*t_ai(b, n_3, l, p)*t_abij(a, n_2 , d, q, k, r, j, n_1)  #69
    term_69 = result - previous_result
    previous_result = result
    print('term69', term_69)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = L_m

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for l in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_n_3)  & set(L_o) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_1)  & set(L_n_3)  & set(L_o)  & set(L_q)  & set(L_p) )
                            for r in domain_r:
                                for k in range(n_occupied):
                                    result += -1*v_pqrs(k, r, l, p, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, i, m)*t_ai(a, n_2, l, p)*t_abij(b, n_3 , d, q, k, r, j, n_1)  #70
    term_70 = result - previous_result
    previous_result = result
    print('term70', term_70)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_n_1

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for l in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_m)  & set(L_o)  & set(L_n_2) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_2)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_m) )
                            for r in domain_r:
                                for k in range(n_occupied):
                                    result += -1*v_pqrs(k, r, l, p, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, j, n_1)*t_ai(b, n_3, l, p)*t_abij(a, n_2 , d, q, k, r, i, m)  #71
    term_71 = result - previous_result
    previous_result = result
    print('term71', term_71)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = L_n_1

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for l in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_3)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_n_3)  & set(L_o)  & set(L_q)  & set(L_p)  & set(L_m) )
                            for r in domain_r:
                                for k in range(n_occupied):
                                    result += 1*v_pqrs(k, r, l, p, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, j, n_1)*t_ai(a, n_2, l, p)*t_abij(b, n_3 , d, q, k, r, i, m)  #72
    term_72 = result - previous_result
    previous_result = result
    print('term72', term_72)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = L_m

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_n_3)  & set(L_q)  & set(L_o) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += 0.25*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, i, m)*t_ai(a, n_2, k, p)*t_ai(d, q, j, n_1)*t_ai(b, n_3, l, r)  #73
    term_73 = result - previous_result
    previous_result = result
    print('term73', term_73)




    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_m

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_n_1)  & set(L_o) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_q)  & set(L_n_2)  & set(L_o) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += -0.25*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, i, m)*t_ai(b, n_3, k, p)*t_ai(d, q, j, n_1)*t_ai(a, n_2, l, r)  #74
    term_74 = result - previous_result
    previous_result = result
    print('term74', term_74)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    domain_o = L_n_1

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_o)  & set(L_n_2) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_n_3)  & set(L_q)  & set(L_o) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += -0.25*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, j, n_1)*t_ai(a, n_2, k, p)*t_ai(d, q, i, m)*t_ai(b, n_3, l, r)  #75
    term_75 = result - previous_result
    previous_result = result
    print('term75', term_75)




    L_n_1 = [i for i in range(n_1-L, n_1+L+1)]
    L_n_3 = [i for i in range(n_3-L, n_3+L+1)]
    L_m = [i for i in range(m-L, m+L+1)]
    L_n_2 = [i for i in range(n_2-L, n_2+L+1)]
    domain_o = L_n_1

    for o in domain_o:
        for c in range(n_virtuals):
            L_o = [i for i in range(o-L, o+L+1)]
            domain_p = list( set(L_n_3)  & set(L_o) )
            for p in domain_p:
                for k in range(n_occupied):
                    L_p = [i for i in range(p-L, p+L+1)]
                    domain_q = list( set(L_p)  & set(L_o)  & set(L_m) )
                    for q in domain_q:
                        for d in range(n_virtuals):
                            L_q = [i for i in range(q-L, q+L+1)]
                            domain_r = list( set(L_p)  & set(L_q)  & set(L_n_2)  & set(L_o) )
                            for r in domain_r:
                                for l in range(n_occupied):
                                    result += 0.25*v_pqrs(k, p, l, r, c + n_occupied, o, d + n_occupied, q)*t_ai(c, o, j, n_1)*t_ai(b, n_3, k, p)*t_ai(d, q, i, m)*t_ai(a, n_2, l, r)  #76
    term_76 = result - previous_result
    previous_result = result
    print('term76', term_76)



    return(result)
