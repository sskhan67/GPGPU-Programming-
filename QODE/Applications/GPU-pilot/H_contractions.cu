/*   (C) Copyright 2018, 2020 Anthony D. Dutoi and Yuhong Liu
 *
 *   This file is part of Qode.
 *
 *   Qode is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Qode is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Qode.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "PyC_types.h"
#include<cuda.h>
#include "cuda_runtime.h"
#include "stdlib.h"
#include "time.h"
#include "stdio.h"
#include "stdlib.h"
#include <sys/time.h>
#include<omp.h>

extern "C"
 {

Double monomer(PyInt n_orb, Double* Rca, Double* Rccaa, Double* h, Double* V)
	{
	Double H = 0;
	//printf("%d N orbibtal: ",n_orb);
	PyInt p=0;
	for (p=0;  p<n_orb;  p++)
		{
		PyInt q=0;
		for (q=0;  q<n_orb;  q++)
			{
			PyInt r=0;
			for ( r=0;  r<n_orb;  r++)
				{
				PyInt s=0;
				for (s=0;  s<n_orb;  s++)
					{
					H += V[((p*n_orb + q)*n_orb + r)*n_orb + s] * Rccaa[((p*n_orb + q)*n_orb + s)*n_orb + r];
					}
				}
			}
		}
	PyInt p1=0;
	for (p1=0;  p1<n_orb;  p1++)
		{
		PyInt q1=0;
		for (q1=0;  q1<n_orb;  q1++)
			{
			H += h[p1*n_orb + q1] * Rca[p1*n_orb + q1];
			}
		}
	return H;
	}



Double dimer_2min2pls(PyInt n_orb1, PyInt n_orb2, Double* Rcc1, Double* Raa2, Double* V)
	{
	Double H = 0;
	PyInt p1=0;
	for (p1=0;  p1<n_orb1;  p1++)
		{
		PyInt q1=0;
		for (q1=0;  q1<n_orb1;  q1++)
			{
			PyInt r2=0;
			for (r2=0;  r2<n_orb2;  r2++)
				{
				PyInt s2=0;
				for ( s2=0;  s2<n_orb2;  s2++)
					{
					H += V[((p1*n_orb1 + q1)*n_orb2 + r2)*n_orb2 + s2] * Rcc1[p1*n_orb1 + q1] * Raa2[s2*n_orb2 + r2];
					}
				}
			}
		}
	return H;
 
	}


Double dimer_1pls1min(PyInt n_orb1, PyInt n_orb2, Double* Ra1, Double* Rcaa1, Double* Rc2, Double* Rcca2, Double* h, Double* V2221, Double* V2111)
	{
	Double H = 0;
	PyInt p2=0;
	for (p2=0;  p2<n_orb2;  p2++)
		{
		PyInt q2=0;
		for ( q2=0;  q2<n_orb2;  q2++)
			{
			PyInt r2=0;
			for ( r2=0;  r2<n_orb2;  r2++)
				{
				PyInt s1=0;
				for ( s1=0;  s1<n_orb1;  s1++)
					{
					H += V2221[((p2*n_orb2 + q2)*n_orb2 + r2)*n_orb1 + s1] * Rcca2[(q2*n_orb2 + p2)*n_orb2 + r2] * Ra1[s1];
					}
				}
			}
		}
	PyInt p22=0;
	for (p22=0;  p22<n_orb2;  p22++)
		{
		PyInt q1=0;
		for ( q1=0;  q1<n_orb1;  q1++)
			{
			PyInt r1=0;
			for (r1=0;  r1<n_orb1;  r1++)
				{
				PyInt s2=0;
				for (s2=0;  s2<n_orb1;  s2++)
					{
					H += V2111[((p22*n_orb1 + q1)*n_orb1 + r1)*n_orb1 + s2] * Rc2[p22] * Rcaa1[(q1*n_orb1 + s2)*n_orb1 + r1];
					}
				}
			}
		}
	H *= 2;
	PyInt p3=0;
	for (p3=0;  p3<n_orb2;  p3++)
		{
		PyInt q3=0;
		for (q3=0;  q3<n_orb1;  q3++)
			{
			H += h[p3*n_orb1 + q3] * Rc2[p3] * Ra1[q3];
			}
		}
	return H;
	}



Double dimer_2pls2min(PyInt n_orb1, PyInt n_orb2, Double* Raa1, Double* Rcc2, Double* V)
	{
	Double H = 0;
	PyInt p2=0;
	for (p2=0;  p2<n_orb2;  p2++)
		{
		PyInt q2=0;
		for (q2=0;  q2<n_orb2;  q2++)
			{
			PyInt r1=0;
			for (r1=0;  r1<n_orb1;  r1++)
				{
				PyInt s1=0;
				for (s1=0;  s1<n_orb1;  s1++)
					{
					H += V[((p2*n_orb2 + q2)*n_orb1 + r1)*n_orb1 + s1] * Rcc2[p2*n_orb2 + q2] * Raa1[s1*n_orb1 + r1];
					}
				}
			}
		}
	return H;
	}


Double dimer_00(PyInt n_orb1, PyInt n_orb2, Double* Rca1, Double* Rca2, Double* V)
	{
	Double H = 0;
	PyInt p1=0;
	for ( p1=0;  p1<n_orb1;  p1++)
		{
		PyInt q2=0;
		for ( q2=0;  q2<n_orb2;  q2++)
			{
			PyInt r1=0;
			for ( r1=0;  r1<n_orb1;  r1++)
				{
				PyInt s2=0;
				for ( s2=0;  s2<n_orb2;  s2++)
					{
					H += V[((p1*n_orb2 + q2)*n_orb1 + r1)*n_orb2 + s2] * Rca1[p1*n_orb1 + r1] * Rca2[q2*n_orb2 + s2];
					}
				}
			}
		}
	return 4*H;
	}





void dimer_1min1pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Rc1,  Double** Rcca1, Double** Ra2, Double** Rcaa2, Double** h, Double** V1, Double** V2)
{
    
	// define  iterator variables
    int n;
    int p1,r1,q1,s2;
    int p11,q2,r2;
    int p12,q22;

    // Data divides into four parts 

    int first_q=n_elem/4;
    int sec_q=first_q*2;
    int third_q=first_q*3;
    
    
    /* First quarter data */

	for(n=0; n<first_q; n++) {
        int index = i[n]*dim[n]+j[n];
        double tmp = 0.0;

    
        for(p1=0; p1<n_orb1[n]; p1++) {
            for(r1=0; r1<n_orb1[n]; r1++) {
                for(q1=0; q1<n_orb1[n]; q1++) {
                    for( s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V1[n][((p1*n_orb1[n] + q1)*n_orb1[n] + r1)*n_orb2[n] + s2] * Rcca1[n][(q1*n_orb1[n] + p1)*n_orb1[n] + r1] * Ra2[n][s2];
                    }
                }
            }
        }

    
        for(p11=0; p11<n_orb1[n]; p11++) {
            for(q2=0; q2<n_orb2[n]; q2++) {
                for(r2=0; r2<n_orb2[n]; r2++) {
                    for(s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V2[n][((p11*n_orb2[n] + q2)*n_orb2[n] + r2)*n_orb2[n] + s2] * Rc1[n][p11] * Rcaa2[n][(q2*n_orb2[n] + s2)*n_orb2[n] + r2];
                    }
                }
            }
        }

        tmp *= 2;

        
        for(p12=0; p12<n_orb1[n]; p12++) {
            for(q22=0; q22<n_orb2[n]; q22++) {
                tmp += h[n][p12*n_orb2[n]+q22] * Rc1[n][p12] * Ra2[n][q22];
            }
        }
        H[n][index] = tmp*sign[n];
    }
    
	/*  2nd Quarter*/
	

	for(n=first_q; n<sec_q; n++) {
        int index = i[n]*dim[n]+j[n];
        double tmp = 0.0;

        
        for(p1=0; p1<n_orb1[n]; p1++) {
            for(r1=0; r1<n_orb1[n]; r1++) {
                for(q1=0; q1<n_orb1[n]; q1++) {
                    for( s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V1[n][((p1*n_orb1[n] + q1)*n_orb1[n] + r1)*n_orb2[n] + s2] * Rcca1[n][(q1*n_orb1[n] + p1)*n_orb1[n] + r1] * Ra2[n][s2];
                    }
                }
            }
        }

        
        for(p11=0; p11<n_orb1[n]; p11++) {
            for(q2=0; q2<n_orb2[n]; q2++) {
                for(r2=0; r2<n_orb2[n]; r2++) {
                    for(s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V2[n][((p11*n_orb2[n] + q2)*n_orb2[n] + r2)*n_orb2[n] + s2] * Rc1[n][p11] * Rcaa2[n][(q2*n_orb2[n] + s2)*n_orb2[n] + r2];
                    }
                }
            }
        }

        tmp *= 2;

        
        for(p12=0; p12<n_orb1[n]; p12++) {
            for(q22=0; q22<n_orb2[n]; q22++) {
                tmp += h[n][p12*n_orb2[n]+q22] * Rc1[n][p12] * Ra2[n][q22];
            }
        }
        H[n][index] = tmp*sign[n];
    }



	/* Third Quarter */

	for(n=sec_q; n<third_q; n++) {
        int index = i[n]*dim[n]+j[n];
        double tmp = 0.0;

        
        for(p1=0; p1<n_orb1[n]; p1++) {
            for(r1=0; r1<n_orb1[n]; r1++) {
                for(q1=0; q1<n_orb1[n]; q1++) {
                    for( s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V1[n][((p1*n_orb1[n] + q1)*n_orb1[n] + r1)*n_orb2[n] + s2] * Rcca1[n][(q1*n_orb1[n] + p1)*n_orb1[n] + r1] * Ra2[n][s2];
                    }
                }
            }
        }

        
        for(p11=0; p11<n_orb1[n]; p11++) {
            for(q2=0; q2<n_orb2[n]; q2++) {
                for(r2=0; r2<n_orb2[n]; r2++) {
                    for(s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V2[n][((p11*n_orb2[n] + q2)*n_orb2[n] + r2)*n_orb2[n] + s2] * Rc1[n][p11] * Rcaa2[n][(q2*n_orb2[n] + s2)*n_orb2[n] + r2];
                    }
                }
            }
        }

        tmp *= 2;

       
        for(p12=0; p12<n_orb1[n]; p12++) {
            for(q22=0; q22<n_orb2[n]; q22++) {
                tmp += h[n][p12*n_orb2[n]+q22] * Rc1[n][p12] * Ra2[n][q22];
            }
        }
        H[n][index] = tmp*sign[n];
    }




	/* Fourth Quarter */

        for(n=third_q; n<n_elem; n++) {
        int index = i[n]*dim[n]+j[n];
        double tmp = 0.0;

        
        for(p1=0; p1<n_orb1[n]; p1++) {
            for(r1=0; r1<n_orb1[n]; r1++) {
                for(q1=0; q1<n_orb1[n]; q1++) {
                    for( s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V1[n][((p1*n_orb1[n] + q1)*n_orb1[n] + r1)*n_orb2[n] + s2] * Rcca1[n][(q1*n_orb1[n] + p1)*n_orb1[n] + r1] * Ra2[n][s2];
                    }
                }
            }
        }

        
        for(p11=0; p11<n_orb1[n]; p11++) {
            for(q2=0; q2<n_orb2[n]; q2++) {
                for(r2=0; r2<n_orb2[n]; r2++) {
                    for(s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V2[n][((p11*n_orb2[n] + q2)*n_orb2[n] + r2)*n_orb2[n] + s2] * Rc1[n][p11] * Rcaa2[n][(q2*n_orb2[n] + s2)*n_orb2[n] + r2];
                    }
                }
            }
        }

        tmp *= 2;

        
        for(p12=0; p12<n_orb1[n]; p12++) {
            for(q22=0; q22<n_orb2[n]; q22++) {
                tmp += h[n][p12*n_orb2[n]+q22] * Rc1[n][p12] * Ra2[n][q22];
            }
        }
        H[n][index] = tmp*sign[n];
    }



    
}


void monomer_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb, Double** Rca,  Double** Rccaa, Double** h, Double** V)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = monomer(n_orb[n], Rca[n], Rccaa[n], h[n], V[n]);
		}
	return;
	}

void dimer_2min2pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Rcc1, Double** Raa2, Double** V)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = dimer_2min2pls(n_orb1[n], n_orb2[n], Rcc1[n], Raa2[n], V[n]);
		}
	return;
	}


void dimer_00_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Rca1, Double** Rca2, Double** V)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = dimer_00(n_orb1[n], n_orb2[n], Rca1[n], Rca2[n], V[n]);
		}
	return;
	}

void dimer_1pls1min_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Ra1,  Double** Rcaa1, Double** Rc2, Double** Rcca2, Double** h, Double** V1, Double** V2)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = sign[n] * dimer_1pls1min(n_orb1[n], n_orb2[n], Ra1[n], Rcaa1[n], Rc2[n], Rcca2[n], h[n], V1[n], V2[n]);
		}
	return;
	}

void dimer_2pls2min_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Raa1, Double** Rcc2, Double** V)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = dimer_2pls2min(n_orb1[n], n_orb2[n], Raa1[n], Rcc2[n], V[n]);
		}
	return;
	}


}
