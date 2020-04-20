
#from tensors import *


#L = 1
def delta_t2(i, j, a, b, n_virtuals, n_occupied, f_pq, v_pqrs, t_ai, t_abij):
    result = 0


    ##print("=================================")
    ##print("i=", i, "j=", j, "a=", a, "b=", b)
    ##print("=================================")
    ##print(t_abij.amplitudes_container, 't2 amp')
    #print(t_abij.amplitudes_container, 'slice of t2')

    result += 1*v_pqrs(a + n_occupied, b + n_occupied, i, j)  #1
    term_1 = result 
    previous_result = result
    #print('term1 =', term_1)

    for c in range(n_virtuals):
            result += 1*f_pq(b + n_occupied, c + n_occupied)*t_abij(a ,c ,i ,j)  #2
    term_2 = result - previous_result
    previous_result = result
    #print('term2=', term_2)

    for c in range(n_virtuals):
            result += -1*f_pq(a + n_occupied, c + n_occupied)*t_abij(b ,c ,i ,j) #3
    term_3 = result - previous_result
    previous_result = result
    #print('term3=', term_3)

    for k in range(n_occupied):
            result += -1*f_pq(k, j)*t_abij(a ,b ,i ,k)  #4
    term_4 = result - previous_result
    previous_result = result
    #print('term4=', term_4)

    for k in range(n_occupied):
            result += 1*f_pq(k, i)*t_abij(a ,b ,j ,k)   #5
    term_5 = result - previous_result
    previous_result = result
    #print('term5=', term_5)

    for k in range(n_occupied):
            for l in range(n_occupied):
                    result += 0.5*v_pqrs(k, l, i, j)*t_abij(a, b, k, l)   #6  
    term_6  = result - previous_result
    previous_result = result
    #print('term6 =', term_6 )

    for d in range(n_virtuals):
            for c in range(n_virtuals):
                    result += 0.5*v_pqrs(a + n_occupied, b + n_occupied, c + n_occupied, d + n_occupied)*t_abij(c, d, i, j)   #7
    term_7 = result - previous_result
    previous_result = result
    #print('term7=', term_7)

    #print('-----------------------------------------')
    for k in range(n_occupied):
            for c in range(n_virtuals):
                    #if v_pqrs(k, b + n_occupied, c + n_occupied, j) != 0.0 and t_abij(a ,c ,i ,k) != 0.0:
                        #print('v=', v_pqrs(k, b + n_occupied, c + n_occupied, j), 't=', t_abij(a ,c ,i ,k), 'multiplied in loop')
                        #print('a=', a, 'c=', c, 'i=', i, 'k=', k)
                    result += 1*v_pqrs(k, b + n_occupied, c + n_occupied, j)*t_abij(a ,c ,i ,k) #8
    term_8 = result - previous_result
    previous_result = result
    #print('term8=', term_8)
    #print('-----------------------------------------')


    for k in range(n_occupied):
            for c in range(n_virtuals):
                    result += -1*v_pqrs(k, b + n_occupied, c + n_occupied, i )*t_abij(a  ,c ,j ,k) #9
    term_9 = result - previous_result
    previous_result = result
    #print('term9=', term_9)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    result += -1*v_pqrs(k, a + n_occupied, c + n_occupied, j)*t_abij(b ,c ,i ,k)   #10
    term_10 = result - previous_result
    previous_result  = result
    #print('term10 =',term_10)



    for k in range(n_occupied):
            for c in range(n_virtuals):
                    result += 1*v_pqrs(k, a + n_occupied, c + n_occupied, i)*t_abij(b ,c ,j ,k)    #11
    term_11= result - previous_result
    previous_result = result
    #print('term11 =', term_11)

    for c in range(n_virtuals):
            result += 1*v_pqrs(a + n_occupied, b + n_occupied, c + n_occupied, j)*t_ai(c, i)    #12
    term_12= result - previous_result
    previous_result = result
    #print('term12 =', term_12)

    for c in range(n_virtuals):
            result += -1*v_pqrs(a + n_occupied, b + n_occupied, c + n_occupied, i)*t_ai(c, j)   #13 mistake
    term_13= result - previous_result
    previous_result = result
    #print('term13 =', term_13)

    for k in range(n_occupied):
            result += -1*v_pqrs(k, b + n_occupied,  i, j)*t_ai(a,k)   #14
    term_14= result - previous_result
    previous_result = result
    #print('term14 =', term_14)

    for k in range(n_occupied):
            result += 1*v_pqrs(k, a + n_occupied, i, j)*t_ai(b, k)   #15 mistake
    term_15= result - previous_result
    previous_result = result
    #print('term15 =', term_15)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += 0.5*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(a ,c ,i ,k)*t_abij(d ,b ,l ,j)   #16
    term_16= result - previous_result
    previous_result = result
    #print('term16 =', term_16)
    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += -0.5*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(b ,c ,i ,k)*t_abij(d ,a ,l ,j)  #17
    term_17= result - previous_result
    previous_result = result
    #print('term17 =', term_17)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += -0.5*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(a ,c ,j ,k)*t_abij(d ,b ,l ,i)  #18
    term_18= result - previous_result
    previous_result = result
    #print('term18 =', term_18)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += 0.5*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(b ,c ,j ,k)*t_abij(d ,a ,l ,i)  #19
    term_19= result - previous_result
    previous_result = result
    #print('term19 =', term_19)

    for d in range(n_virtuals):
            for c in range(n_virtuals):
                    for k in range(n_occupied):
                            for l in range(n_occupied):
                                    result += 0.25*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(c ,d ,i ,j)*t_abij(a ,b ,k ,l)  #20
    term_20= result - previous_result
    previous_result = result
    #print('term20 =', term_20)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += -0.5*v_pqrs(k, l,  c + n_occupied, d + n_occupied)*t_abij(a ,c ,i ,j)*t_abij(b ,d ,k ,l)  #21
    term_21= result - previous_result
    previous_result = result
    #print('term21 =', term_21)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += 0.5*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(b ,c ,i ,j)*t_abij(a ,d ,k ,l)  #22
    term_22= result - previous_result
    previous_result = result
    #print('term22 =', term_22)


    for k in range(n_occupied):
            for d in range(n_virtuals):
                    for c in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += -0.5*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(a ,b ,i ,k)*t_abij(c ,d ,j ,l)  #23
    term_23= result - previous_result
    previous_result = result
    #print('term23 =', term_23)

    for k in range(n_occupied):
            for d in range(n_virtuals):
                    for c in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += 0.5*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(a ,b ,j ,k)*t_abij(c ,d ,i ,l)  #24
    term_24= result - previous_result
    previous_result = result
    #print('term24 =', term_24)

    for k in range(n_occupied):
            for l in range(n_occupied):
                    result += 0.5*v_pqrs(k, l, i, j)*t_ai(a, k)*t_ai(b, l)  #25
    term_25= result - previous_result
    previous_result = result
    #print('term25 =', term_25)


    for k in range(n_occupied):
            for l in range(n_occupied):
                    result += -0.5*v_pqrs(k, l, i, j)*t_ai(b, k)*t_ai(a, l)  #26
    term_26= result - previous_result
    previous_result = result
    #print('term26 =', term_26)


    for c in range(n_virtuals):
            for d in range(n_virtuals):
                    result += 0.5*v_pqrs(a + n_occupied, b + n_occupied, c + n_occupied, d + n_occupied)*t_ai(c, i)*t_ai(d, j)   #27
    term_27= result - previous_result
    previous_result = result
    #print('term27 =', term_27)


    for c in range(n_virtuals):
            for d in range(n_virtuals):
                    result += -0.5*v_pqrs(a + n_occupied, b + n_occupied, c + n_occupied, d + n_occupied)*t_ai(c, j)*t_ai(d, i)   #28
    term_28= result - previous_result
    previous_result = result
    #print('term28 =', term_28)


    for k in range(n_occupied):
            for c in range(n_virtuals):
                    result += -1*v_pqrs(k, b + n_occupied, i, c + n_occupied)*t_ai(a, k)*t_ai(c, j)  #29
    term_29= result - previous_result
    previous_result = result
    #print('term29 =', term_29)


    for k in range(n_occupied):
            for c in range(n_virtuals):
                    result += 1*v_pqrs(k, a + n_occupied, i, c + n_occupied)*t_ai(b, k)*t_ai(c, j)   #30
    term_30= result - previous_result
    previous_result = result
    #print('term30 =', term_30)


    for k in range(n_occupied):
            for c in range(n_virtuals):
                    result += 1*v_pqrs(k, b + n_occupied, j, c + n_occupied)*t_ai(a, k)*t_ai(c, i)   #31
    term_31= result - previous_result
    previous_result = result
    #print('term31 =', term_31)


    for k in range(n_occupied):
            for c in range(n_virtuals):
                    result += -1*v_pqrs(k, a + n_occupied, j, c + n_occupied)*t_ai(b, k)*t_ai(c, i)  #32
    term_32= result - previous_result
    previous_result = result
    #print('term32 =', term_32)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    result += 1*f_pq(k, c + n_occupied)*t_ai(a, k)*t_abij(b ,c ,i ,j)     #33
    term_33= result - previous_result
    previous_result = result
    #print('term33 =', term_33)


    for k in range(n_occupied):
            for c in range(n_virtuals):
                    result += -1*f_pq(k, c + n_occupied)*t_ai(b, k )*t_abij(a ,c ,i ,j)   #34
    term_34= result - previous_result
    previous_result = result
    #print('term34 =', term_34)


    for c in range(n_virtuals):
            for k in range(n_occupied):
                    result += 1*f_pq(k, c + n_occupied)*t_ai(c ,i)*t_abij(a ,b ,j ,k)    #35
    term_35= result - previous_result
    previous_result = result
    #print('term35 =', term_35)


    for c in range(n_virtuals):
            for k in range(n_occupied):
                    result += -1*f_pq(k, c + n_occupied)*t_ai(c ,j)*t_abij(a ,b ,i ,k)   #36
    term_36= result - previous_result
    previous_result = result
    #print('term36 =', term_36)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for l in range(n_occupied):
                            result += -1*v_pqrs(k, l, c + n_occupied, i)*t_ai(c ,k )*t_abij(a ,b ,l ,j)   #37
    term_37= result - previous_result
    previous_result = result
    #print('term37 =', term_37)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for l in range(n_occupied):
                            result += 1*v_pqrs(k, l, c + n_occupied, j)*t_ai(c ,k )*t_abij(a ,b ,l ,i)    #38
    term_38= result - previous_result
    previous_result = result
    #print('term38 =', term_38)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for d in range(n_virtuals):
                            result += 1*v_pqrs(k, a + n_occupied, c + n_occupied, d + n_occupied)*t_ai(c ,k )*t_abij(d  ,b ,i ,j)  #39
    term_39= result - previous_result
    previous_result = result
    #print('term39 =', term_39)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for d in range(n_virtuals):
                            result += -1*v_pqrs(k, b + n_occupied, c + n_occupied, d + n_occupied)*t_ai(c ,k )*t_abij(d  ,a ,i ,j)  #40
    term_40= result - previous_result
    previous_result = result
    #print('term40 =', term_40)


    for d in range(n_virtuals):
            for k in range(n_occupied):
                    for c in range(n_virtuals):
                            result += 1*v_pqrs(a + n_occupied, k, d + n_occupied, c + n_occupied)*t_ai(d ,i)*t_abij(b  ,c ,j ,k)   #41
    term_41= result - previous_result
    previous_result = result
    #print('term41 =', term_41)


    for d in range(n_virtuals):
            for k in range(n_occupied):
                    for c in range(n_virtuals):
                            result += -1*v_pqrs(b + n_occupied, k, d + n_occupied, c + n_occupied)*t_ai(d ,i)*t_abij(a  ,c ,j ,k)  #42
    term_42= result - previous_result
    previous_result = result
    #print('term42 =', term_42)


    for d in range(n_virtuals):
            for k in range(n_occupied):
                    for c in range(n_virtuals):
                            result += -1*v_pqrs(a + n_occupied, k, d + n_occupied, c + n_occupied)*t_ai(d ,j)*t_abij(b  ,c ,i ,k)  #43
    term_43= result - previous_result
    previous_result = result
    #print('term43 =', term_43)


    for d in range(n_virtuals):
            for k in range(n_occupied):
                    for c in range(n_virtuals):
                            result += 1*v_pqrs(b + n_occupied, k, d + n_occupied, c + n_occupied)*t_ai(d ,j)*t_abij(a  ,c ,i ,k)   #44
    term_44= result - previous_result
    previous_result = result
    #print('term44 =', term_44)


    for l in range(n_occupied):
            for k in range(n_occupied):
                    for c in range(n_virtuals):
                            result += 1*v_pqrs(k, l, i, c + n_occupied)*t_ai(a, l )*t_abij(b  ,c ,j ,k)  #45
    term_45= result - previous_result
    previous_result = result
    #print('term45 =', term_45)


    for l in range(n_occupied):
            for k in range(n_occupied):
                    for c in range(n_virtuals):
                            result += -1*v_pqrs(k, l, i, c + n_occupied)*t_ai(b, l )*t_abij(a  ,c ,j ,k)  #46
    term_46= result - previous_result
    previous_result = result
    #print('term46 =', term_46)


    for l in range(n_occupied):
            for k in range(n_occupied):
                    for c in range(n_virtuals):
                            result += -1*v_pqrs(k, l, j, c + n_occupied)*t_ai(a, l )*t_abij(b  ,c ,i ,k)  #47
    term_47= result - previous_result
    previous_result = result
    #print('term47 =', term_47)


    for l in range(n_occupied):
            for k in range(n_occupied):
                    for c in range(n_virtuals):
                            result += 1*v_pqrs(k, l, j, c + n_occupied)*t_ai(b, l )*t_abij(a  ,c ,i ,k)   #48
    term_48= result - previous_result
    previous_result = result
    #print('term48 =', term_48)


    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for l in range(n_occupied):
                            result += 0.5*v_pqrs(k, l, c + n_occupied, j)*t_ai(c ,i)*t_abij(a  ,b ,k ,l)  #49
    term_49= result - previous_result
    previous_result = result
    #print('term49 =', term_49)


    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for l in range(n_occupied):
                            result += -0.5*v_pqrs(k, l, c + n_occupied, i)*t_ai(c ,j)*t_abij(a  ,b ,k ,l)  #50
    term_50= result - previous_result
    previous_result = result
    #print('term50 =', term_50)


    for k in range(n_occupied):
            for d in range(n_virtuals):
                    for c in range(n_virtuals):
                            result += -0.5*v_pqrs(k, b + n_occupied, c + n_occupied, d + n_occupied)*t_ai(a, k )*t_abij(c  ,d ,i ,j)  #51
    term_51= result - previous_result
    previous_result = result
    #print('term51 =', term_51)


    for k in range(n_occupied):
            for d in range(n_virtuals):
                    for c in range(n_virtuals):
                            result += 0.5*v_pqrs(k, a + n_occupied, c + n_occupied, d + n_occupied)*t_ai(b, k )*t_abij(c  ,d ,i ,j)   #52
    term_52= result - previous_result
    previous_result = result
    #print('term52 =', term_52)


    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            result += -0.5*v_pqrs(k, b + n_occupied, c + n_occupied, d + n_occupied)*t_ai(c ,i)*t_ai(a ,k)*t_ai(d ,j)  #53
    term_53= result - previous_result
    previous_result = result
    #print('term53 =', term_53)


    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            result += 0.5*v_pqrs(k, a + n_occupied, c + n_occupied, d + n_occupied)*t_ai(c ,i)*t_ai(b ,k)*t_ai(d ,j)  #54
    term_54= result - previous_result
    previous_result = result
    #print('term54 =', term_54)
    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            result += 0.5*v_pqrs(k, b + n_occupied, c + n_occupied, d + n_occupied)*t_ai(c ,j)*t_ai(a ,k)*t_ai(d ,i)  #55 mistake
    term_55= result - previous_result
    previous_result = result
    #print('term55 =', term_55)


    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            result += -0.5*v_pqrs(k, a + n_occupied, c + n_occupied, d + n_occupied)*t_ai(c ,j)*t_ai(b ,k)*t_ai(d ,i)  #56
    term_56= result - previous_result
    previous_result = result
    #print('term56 =', term_56)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for l in range(n_occupied):
                            result += 0.5*v_pqrs(k, l, c + n_occupied, j)*t_ai(c ,i)*t_ai(a ,k)*t_ai(b ,l)   #57
    term_57= result - previous_result
    previous_result = result
    #print('term57 =', term_57)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for l in range(n_occupied):
                            result += -0.5*v_pqrs(k, l, c + n_occupied, j)*t_ai(c ,i)*t_ai(b ,k)*t_ai(a ,l)  #58  here
    term_re= result - previous_result
    previous_result = result
    #print('termre =', term_re)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for l in range(n_occupied):
                            result += -0.5*v_pqrs(k, l, c + n_occupied, i)*t_ai(c ,j)*t_ai(a ,k)*t_ai(b ,l)  #59
    term_59= result - previous_result
    previous_result = result
    #print('term59 =', term_59)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for l in range(n_occupied):
                            result += 0.5*v_pqrs(k, l, c + n_occupied, i)*t_ai(c ,j)*t_ai(b ,k)*t_ai(a ,l)  #60
    term_60= result - previous_result
    previous_result = result
    #print('term60 =', term_60)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += -1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,k )*t_ai(d ,i)*t_abij(a ,b ,l ,j)   #61
    term_61= result - previous_result
    previous_result = result
    #print('term61 =', term_61)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += 1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,k )*t_ai(d ,j)*t_abij(a ,b ,l ,i)   #62
    term_62= result - previous_result
    previous_result = result
    #print('term62 =', term_62)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for l in range(n_occupied):
                            for d in range(n_virtuals):
                                    result += -1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,k )*t_ai(a, l)*t_abij(d  ,b ,i ,j)  #63
    term_63= result - previous_result
    previous_result = result
    #print('term63 =', term_63)

    for k in range(n_occupied):
            for c in range(n_virtuals):
                    for l in range(n_occupied):
                            for d in range(n_virtuals):
                                    result += 1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,k )*t_ai(b, l)*t_abij(d  ,a ,i ,j)  #64
    term_64= result - previous_result
    previous_result = result
    #print('term64 =', term_64)

    for c in range(n_virtuals):
            for d in range(n_virtuals):
                    for k in range(n_occupied):
                            for l in range(n_occupied):
                                    result += 0.25*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,i)*t_ai(d ,j)*t_abij(a  ,b ,k ,l )  #65
    term_65= result - previous_result
    previous_result = result
    #print('term65 =', term_65)

    for c in range(n_virtuals):
            for d in range(n_virtuals):
                    for k in range(n_occupied):
                            for l in range(n_occupied):
                                    result += -0.25*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,j)*t_ai(d ,i)*t_abij(a  ,b ,k ,l )  #66
    term_66= result - previous_result
    previous_result = result
    #print('term66 =', term_66)

    for k in range(n_occupied):
            for l in range(n_occupied):
                    for d in range(n_virtuals):
                            for c in range(n_virtuals):
                                    result += 0.25*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(a, k )*t_ai(b, l)*t_abij(c  ,d ,i ,j)  #67
    term_67= result - previous_result
    previous_result = result
    #print('term67 =', term_67)

    for k in range(n_occupied):
            for l in range(n_occupied):
                    for d in range(n_virtuals):
                            for c in range(n_virtuals):
                                    result += -0.25*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(b, k )*t_ai(a, l)*t_abij(c  ,d ,i ,j)  #68
    term_68= result - previous_result
    previous_result = result
    #print('term68 =', term_68)

    for c in range(n_virtuals):
            for l in range(n_occupied):
                    for k in range(n_occupied):
                            for d in range(n_virtuals):
                                    result += 1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,i)*t_ai(b, l)*t_abij(a  ,d ,k ,j)   #69
    term_69= result - previous_result
    previous_result = result
    #print('term69 =', term_69)

    for c in range(n_virtuals):
            for l in range(n_occupied):
                    for k in range(n_occupied):
                            for d in range(n_virtuals):
                                    result += -1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,i)*t_ai(a, l)*t_abij(b  ,d ,k ,j)  #70
    term_70= result - previous_result
    previous_result = result
    #print('term70 =', term_70)

    for c in range(n_virtuals):
            for l in range(n_occupied):
                    for k in range(n_occupied):
                            for d in range(n_virtuals):
                                    result += -1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,j)*t_ai(b, l)*t_abij(a  ,d ,k ,i)  #71
    term_71= result - previous_result
    previous_result = result
    #print('term71 =', term_71)

    for c in range(n_virtuals):
            for l in range(n_occupied):
                    for k in range(n_occupied):
                            for d in range(n_virtuals):
                                    result += 1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,j)*t_ai(a, l)*t_abij(b  ,d ,k ,i)   #72
    term_72= result - previous_result
    previous_result = result
    #print('term72 =', term_72)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += 0.25*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,i)*t_ai(a, k)*t_ai(d ,j)*t_ai(b ,l )  #73
    term_73= result - previous_result
    previous_result = result
    #print('term73 =', term_73)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += -0.25*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,i)*t_ai(b, k)*t_ai(d ,j)*t_ai(a ,l )  #74
    term_74= result - previous_result
    previous_result = result
    #print('term74 =', term_74)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += -0.25*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,j)*t_ai(a, k)*t_ai(d ,i)*t_ai(b ,l )  #75
    term_75= result - previous_result
    previous_result = result
    #print('term75 =', term_75)

    for c in range(n_virtuals):
            for k in range(n_occupied):
                    for d in range(n_virtuals):
                            for l in range(n_occupied):
                                    result += 0.25*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,j)*t_ai(b, k)*t_ai(d ,i)*t_ai(a ,l )  #76 mistake
    term_76= result - previous_result
    previous_result = result
    #print('term76 =', term_76)

    ###print(result, 'result final')
    return(result)
