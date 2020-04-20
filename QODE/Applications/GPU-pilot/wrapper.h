#include <cuda.h>


extern void preMalloc(int n_orb1, int n_orb2);
extern double test_wrapper( int n_orb1, int n_orb2, double* Rc1, double* Rcca1, double* Ra2,double* Rcaa2, double* h, double* V1112, double* V1222, int freevariables);

//extern void test_wrapper(float *x, float *y, float *io, int N);
