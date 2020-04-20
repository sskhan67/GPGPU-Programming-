import numpy as np
#from energy_from_amplitudes import H_spatial_to_spin_orb

#from tensors import *


#L = 1
def delta_E(n_virtuals, n_occupied, f_pq, v_pqrs, t_ai, t_abij):
    result = 0


    for a in range(n_virtuals):
        for i in range(n_occupied):
            result += 1*f_pq(i, a + n_occupied)*t_ai(a, i)

    #print("============= CC_MP2 ==============")

    for a in range(n_virtuals):
        for j in range(n_occupied):
            for i in range(n_occupied):
                for b in range(n_virtuals):
                    #if v_pqrs(i, j, a + n_occupied, b + n_occupied)*t_abij(a, b, i ,j) != 0:
                    #    print("i=",i, "j=",j, "a=",a,"b=",b, "V=", v_pqrs(i, j, a + n_occupied, b + n_occupied), "t=", t_abij(a, b, i ,j) )
                    result += 0.25*v_pqrs(i, j, a + n_occupied, b + n_occupied)*t_abij(a, b, i ,j)




 
    for a in range(n_virtuals):
        for i in range(n_occupied):
            for j in range(n_occupied):
                for b in range(n_virtuals):
                     result += 0.5*v_pqrs(i, j, a + n_occupied, b + n_occupied)*t_ai(a, i )*t_ai(b, j)


    return(result)
