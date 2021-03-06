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
#include"wrapper.h"
#include "time.h"
#include "stdio.h"
#include <sys/time.h>

#define PRINT_TIMES 1
#define PRINT_INSIDE_TIMES 0



Double monomer(PyInt n_orb, Double* Rca, Double* Rccaa, Double* h, Double* V)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
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
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_INSIDE_TIMES) printf("monomer finished in %f ms\n", t);
    return H;
    }



Double dimer_2min2pls(PyInt n_orb1, PyInt n_orb2, Double* Rcc1, Double* Raa2, Double* V)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
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
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_INSIDE_TIMES) printf("dimer_2min2pls finished in %f ms\n", t);
    return H;
 
}

static count=0;
Double dimer_1min1pls(PyInt n_orb1, PyInt n_orb2, Double* Rc1, Double* Rcca1, Double* Ra2, Double* Rcaa2, Double* h, Double* V1112, Double* V1222, int freevariables)
	{
	test_wrapper(n_orb1, n_orb2,  Rc1, Rcca1,  Ra2, Rcaa2,  h,  V1112,  V1222,freevariables);
//	count++;
//if(count>10)
	//exit(0);

	}
	
	
/*

Double dimer_1min1pls(PyInt n_orb1, PyInt n_orb2, Double* Rc1, Double* Rcca1, Double* Ra2, Double* Rcaa2, Double* h, Double* V1112, Double* V1222)
    {
    Double H = 0;
    //printf("%d ord1: ",n_orb1);
    //printf(" n_orb2:%d\n ",n_orb2);
    //printf(" Rc1:%f \n",*Rc1);
    //printf(" Rcca1:%f\n ",*Rcca1);
    //printf(" Ra2: %f\n ",*Ra2);
    //printf(" h:%f",*h);
    //printf(" V1112:%f",* V1112);
    PyInt p1=0;
    for ( p1=0;  p1<n_orb1;  p1++)
        {
        PyInt q1=0;
        for ( q1=0;  q1<n_orb1;  q1++)
            {
            PyInt r1=0;
            for ( r1=0;  r1<n_orb1;  r1++)
                {
                PyInt s2=0;
                for ( s2=0;  s2<n_orb2;  s2++)
                    {
                    H += V1112[((p1*n_orb1 + q1)*n_orb1 + r1)*n_orb2 + s2] * Rcca1[(q1*n_orb1 + p1)*n_orb1 + r1] * Ra2[s2];
                    }
                }
            }
        }
    PyInt p11=0;
    for (p11=0;  p11<n_orb1;  p11++)
        {
        PyInt q2=0;
        for (q2=0;  q2<n_orb2;  q2++)
            {
            PyInt r2=0;
            for (r2=0;  r2<n_orb2;  r2++)
                {
                PyInt s2=0;
                for ( s2=0;  s2<n_orb2;  s2++)
                    {
                    H += V1222[((p11*n_orb2 + q2)*n_orb2 + r2)*n_orb2 + s2] * Rc1[p1] * Rcaa2[(q2*n_orb2 + s2)*n_orb2 + r2];
                    }
                }
            }
        }
    H *= 2;
    PyInt p12=0;
    for (p12=0;  p12<n_orb1;  p12++)
        {
        PyInt q22=0;
        for (q22=0;  q22<n_orb2;  q22++)
            {
            H += h[p12*n_orb2 + q22] * Rc1[p12] * Ra2[q22];
            }
        }
    //printf("%f H: ",H);
    return H;
    }

*/
Double dimer_00(PyInt n_orb1, PyInt n_orb2, Double* Rca1, Double* Rca2, Double* V)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
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
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_INSIDE_TIMES) printf("dimer_00 finished in %f ms\n", t);
    return 4*H;
    }



Double dimer_1pls1min(PyInt n_orb1, PyInt n_orb2, Double* Ra1, Double* Rcaa1, Double* Rc2, Double* Rcca2, Double* h, Double* V2221, Double* V2111)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
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
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_INSIDE_TIMES) printf("dimer_1pls1min finished in %f ms\n", t);
    return H;
    }



Double dimer_2pls2min(PyInt n_orb1, PyInt n_orb2, Double* Raa1, Double* Rcc2, Double* V)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
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
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_INSIDE_TIMES) printf("dimer_2pls2min finished in %f ms\n", t);
    return H;
    }





void monomer_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb, Double** Rca,  Double** Rccaa, Double** h, Double** V)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
    PyInt n=0;
    for ( n=0;  n<n_elem;  n++)
        {
        PyInt index = i[n]*dim[n] + j[n];
        H[n][index] = monomer(n_orb[n], Rca[n], Rccaa[n], h[n], V[n]);
        }
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_TIMES) printf("monomer_loop finished in %f ms\n", t);
    }

void dimer_2min2pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Rcc1, Double** Raa2, Double** V)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
    PyInt n=0;
    for ( n=0;  n<n_elem;  n++)
        {
        PyInt index = i[n]*dim[n] + j[n];
        H[n][index] = dimer_2min2pls(n_orb1[n], n_orb2[n], Rcc1[n], Raa2[n], V[n]);
        }
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_TIMES) printf("dimer_2min2pls_loop finished in %f ms\n", t);
    return;
    }

void dimer_1min1pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Rc1,  Double** Rcca1, Double** Ra2, Double** Rcaa2, Double** h, Double** V1, Double** V2)
	{
    struct timeval start,stop; 
    gettimeofday(&start,0);

	int freevariables=0; // not execute dynamic memory free code block , return H
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = sign[n] * dimer_1min1pls(n_orb1[n], n_orb2[n], Rc1[n], Rcca1[n], Ra2[n], Rcaa2[n], h[n], V1[n], V2[n],freevariables);
		}
 	freevariables=1; // execute dynamic memory allocation , retun 0
	
	dimer_1min1pls(n_orb1[0], n_orb2[0], Rc1[0], Rcca1[0], Ra2[0], Rcaa2[0], h[0], V1[0], V2[0],freevariables);
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_TIMES) printf("dimer_1min1pls_loop finished in %f ms\n", t);
	return;
	}

void dimer_00_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Rca1, Double** Rca2, Double** V)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
    PyInt n=0;
    for ( n=0;  n<n_elem;  n++)
        {
        PyInt index = i[n]*dim[n] + j[n];
        H[n][index] = dimer_00(n_orb1[n], n_orb2[n], Rca1[n], Rca2[n], V[n]);
        }
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_TIMES) printf("dimer_00_loop finished in %f ms\n", t);
    return;
    }

void dimer_1pls1min_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Ra1,  Double** Rcaa1, Double** Rc2, Double** Rcca2, Double** h, Double** V1, Double** V2)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
    PyInt n=0;
    for ( n=0;  n<n_elem;  n++)
        {
        PyInt index = i[n]*dim[n] + j[n];
        H[n][index] = sign[n] * dimer_1pls1min(n_orb1[n], n_orb2[n], Ra1[n], Rcaa1[n], Rc2[n], Rcca2[n], h[n], V1[n], V2[n]);
        }
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_TIMES) printf("dimer_1pls1min_loop finished in %f ms\n", t);
    return;
    }

void dimer_2pls2min_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Raa1, Double** Rcc2, Double** V)
    {
    struct timeval start,stop; 
    gettimeofday(&start,0);
    PyInt n=0;
    for ( n=0;  n<n_elem;  n++)
        {
        PyInt index = i[n]*dim[n] + j[n];
        H[n][index] = dimer_2pls2min(n_orb1[n], n_orb2[n], Raa1[n], Rcc2[n], V[n]);
        }
    gettimeofday(&stop,0);
    double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
    if(PRINT_TIMES) printf("dimer_2pls2min_loop finished in %f ms\n", t);
    return;
    }



