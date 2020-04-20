
#from tensors import *


#L = 1
def delta_t(a, i, n_virtuals, n_occupied, f_pq, v_pqrs, t_ai, t_abij):
    result = 0

    previous_result = 0
    result += 1*f_pq(a + n_occupied, i) #1
    previous_result = result - previous_result
    #print('term1 =', result)


    for c in range(n_virtuals):
        result += 1*f_pq(a + n_occupied, c + n_occupied)*t_ai(c, i) #2


    for k in range(n_occupied):
        result += -1*f_pq(k, i)*t_ai(a, k) #3
    previous_result = result - previous_result
    #print('term3 =', result)


    for c in range(n_virtuals):
        for k in range(n_occupied):
            result += 1*v_pqrs(k, a + n_occupied, c + n_occupied, i)*t_ai(c, k) #4
    previous_result = result - previous_result
    #print('term4 =', result)



    for c in range(n_virtuals):
        for k in range(n_occupied):
            result += 1*f_pq(k, c + n_occupied)*t_abij(a  ,c ,i,k)  #5
    previous_result = result - previous_result
    #print('term5 =', result)



    for d in range(n_virtuals):
        for c in range(n_virtuals):
            for k in range(n_occupied):
                result += 0.5*v_pqrs(k, a + n_occupied, c + n_occupied, d + n_occupied)*t_abij(c, d, k, i)  #6
    previous_result = result - previous_result
    #print('term6 =', result)




    for c in range(n_virtuals):
        for l in range(n_occupied):
            for k in range(n_occupied):
                result += -0.5*v_pqrs(k, l, c + n_occupied, i)*t_abij(c  ,a ,k ,l)  #7
    previous_result = result - previous_result
    #print('term7 =', result)



    for c in range(n_virtuals):
        for k in range(n_occupied):
            result += -1*f_pq(k, c + n_occupied)*t_ai(c, i )*t_ai(a, k)  #8
    previous_result = result - previous_result
    #print('term8 =', result)


    for c in range(n_virtuals):
        for k in range(n_occupied):
            for l in range(n_occupied):
                result += -1*v_pqrs(k, l, c + n_occupied, i)*t_ai(c ,k)*t_ai(a,l)  #9
    previous_result = result - previous_result
    #print('term9 =', result)



    for c in range(n_virtuals):
        for k in range(n_occupied):
            for d in range(n_virtuals):
                result += 1*v_pqrs(k, a + n_occupied, c + n_occupied, d + n_occupied)*t_ai(c ,k)*t_ai(d,i)  #10
    previous_result = result - previous_result
    #print('term10 =', result)




    for c in range(n_virtuals):
        for k in range(n_occupied):
            for d in range(n_virtuals):
                for l in range(n_occupied):
                    result += -1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c,k)*t_ai(d ,i)*t_ai(a,l)  #11
    previous_result = result - previous_result
    #print('term11 =', result)




    for c in range(n_virtuals):
        for k in range(n_occupied):
            for d in range(n_virtuals):
                for l in range(n_occupied):
                    result += 1*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_ai(c ,k)*t_abij(d, a, l, i) #12
    previous_result = result - previous_result
    #print('term12 =', result)



    for d in range(n_virtuals):
        for c in range(n_virtuals):
            for k in range(n_occupied):
                for l in range(n_occupied):
                    result += -0.5*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(c  ,d ,k ,i)*t_ai(a,l)  #13
    previous_result = result - previous_result
    #print('term13 =', result)






    for c in range(n_virtuals):
        for l in range(n_occupied):
            for k in range(n_occupied):
                for d in range(n_virtuals):
                    result += -0.5*v_pqrs(k, l, c + n_occupied, d + n_occupied)*t_abij(c  ,a ,k ,l)*t_ai(d ,i) #14
    previous_result = result - previous_result
    #print('term14 =', result)
    #print('result=', result) 
    return result
