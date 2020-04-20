from tensors import *
from generated_chain_code import delta_t
#from generated_chain_code import delta_E

a = 0
n = 0
i = 0
m = 0
L = limit

print('The test result is', delta_t(a, n, i, m, L, n_virtuals, n_occupied, f_pq, v_pqrs, t_ai, t_abij))
print('The equations ran succesfully')
