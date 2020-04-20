import numpy as np
import matplotlib
import matplotlib.pyplot as plt

value_amplitude = 0.77232514

data = np.load('../integrals_files/data.npy')
x_domain = data[0]
term1    = data[1]
term2    = data[2]
term3    = data[3]
term4    = data[4]
term5    = data[5]
term6    = data[6]
term7    = data[7]
term8    = data[8]
term9    = data[9]
term10   = data[10]
term11   = data[11]
term12   = data[12]
term13   = data[13]
term14   = data[14]
total = data[-1]

plt.xlabel('Amplitude value introduced')
plt.ylabel('Equation result')
plt.title('Eq result vs value introduced', fontsize=14, fontweight='bold')
plt.axvline(x=value_amplitude, label='psi4 amp = {}'.format(value_amplitude))
plt.plot(x_domain, total,  label='TOTAL SUM', marker='o', color='r')
plt.plot(x_domain, term1,  label='term 1',  linestyle=':')
plt.plot(x_domain, term2,  label='term 2',  linestyle='-')
plt.plot(x_domain, term3,  label='term 3',  linestyle='-.')
plt.plot(x_domain, term4,  label='term 4',  linestyle='-')
plt.plot(x_domain, term5,  label='term 5',  marker='.')
plt.plot(x_domain, term6,  label='term 6',  linestyle='--')
plt.plot(x_domain, term7,  label='term 7',  linestyle='-', marker='x')
plt.plot(x_domain, term8,  label='term 8',  linestyle='--', color='c')
plt.plot(x_domain, term9,  label='term 9',  linestyle='-')
plt.plot(x_domain, term10, label='term 10', linestyle='-')
plt.plot(x_domain, term11, label='term 11', linestyle=':')
plt.plot(x_domain, term12, label='term 12', linestyle='-.')
plt.plot(x_domain, term13, label='term 13', marker='.')
plt.plot(x_domain, term14, label='term 14', marker=',')
#plt.plot(x_domain, total, label='term1','ro', x_domain, term1, 'b--', x_domain, term2, 'g--', x_domain, term3, 'c--', x_domain, term4, 'm--', x_domain, term5, 'y--', x_domain, term6, 'k--', x_domain, term7, 'w--', x_domain, term8, 'b:', x_domain, term9, 'g:', x_domain, term10, 'c:', x_domain, term11, 'm:', x_domain, term12, 'y:', x_domain, term13, 'k:', x_domain, term14, 'w:')
plt.legend(bbox_to_anchor=(1, 1), loc='best', borderaxespad=0. )
plt.tight_layout(pad=7)
plt.show()
