/*   (C) Copyright 2018, 2019 Anthony D. Dutoi and Yuhong Liu
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
#include "sys/time.h"
#include "stdio.h"
#include "stdlib.h"
#include "cuda.h"
#include "wrapper.h"
#define TIME_INSIDE_FUNCTION 0
#define TIME_LOOP 1
double total_time;

void cutest_f()
{
	int N = 1 << 20;
    float *x = 		malloc(sizeof(float) *N);
    float *y = 		malloc(sizeof(float) *N);
	float *out =	malloc(sizeof(float) *N);
	test_wrapper(x,y,out,N); // declared in header file; defined in .cu
}
Double monomer(PyInt n_orb, Double* Rca, Double* Rccaa, Double* h, Double* V)
{
	struct timeval t1, t2;
	gettimeofday(&t1, NULL);
	Double H = 0;
	for (PyInt p=0;  p<n_orb;  p++)
		{
		for (PyInt q=0;  q<n_orb;  q++)
			{
			for (PyInt r=0;  r<n_orb;  r++)
				{
				for (PyInt s=0;  s<n_orb;  s++)
					{
					H += V[((p*n_orb + q)*n_orb + r)*n_orb + s] * Rccaa[((p*n_orb + q)*n_orb + s)*n_orb + r];
					}
				}
			}
		}
	for (PyInt p=0;  p<n_orb;  p++)
		{
		for (PyInt q=0;  q<n_orb;  q++)
			{
			H += h[p*n_orb + q] * Rca[p*n_orb + q];
			}
		}
	if(TIME_INSIDE_FUNCTION) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		printf("monomer finished execution in %f seconds\n", t);
	}
	return H;
	}



Double dimer_2min2pls(PyInt n_orb1, PyInt n_orb2, Double* Rcc1, Double* Raa2, Double* V)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	Double H = 0;
	for (PyInt p1=0;  p1<n_orb1;  p1++)
	{
		for (PyInt q1=0;  q1<n_orb1;  q1++)
		{
			for (PyInt r2=0;  r2<n_orb2;  r2++)
			{
				for (PyInt s2=0;  s2<n_orb2;  s2++)
				{
					H += V[((p1*n_orb1 + q1)*n_orb2 + r2)*n_orb2 + s2] * Rcc1[p1*n_orb1 + q1] * Raa2[s2*n_orb2 + r2];
				}
			}
		}
	}
	if(TIME_INSIDE_FUNCTION) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		printf("dimer_2min2pls finished execution in %f seconds\n", t);
	}
	return H;
}



Double dimer_1min1pls(PyInt n_orb1, PyInt n_orb2, Double* Rc1, Double* Rcca1, Double* Ra2, Double* Rcaa2, Double* h, Double* V1112, Double* V1222)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	Double H = 0;
	for (PyInt p1=0;  p1<n_orb1;  p1++)
	{
		for (PyInt q1=0;  q1<n_orb1;  q1++)
		{
			for (PyInt r1=0;  r1<n_orb1;  r1++)
			{
				for (PyInt s2=0;  s2<n_orb2;  s2++)
				{
					H += V1112[((p1*n_orb1 + q1)*n_orb1 + r1)*n_orb2 + s2] * Rcca1[(q1*n_orb1 + p1)*n_orb1 + r1] * Ra2[s2];
				}
			}
		}
	}
	for (PyInt p1=0;  p1<n_orb1;  p1++)
	{
		for (PyInt q2=0;  q2<n_orb2;  q2++)
		{
			for (PyInt r2=0;  r2<n_orb2;  r2++)
			{
				for (PyInt s2=0;  s2<n_orb2;  s2++)
				{
					H += V1222[((p1*n_orb2 + q2)*n_orb2 + r2)*n_orb2 + s2] * Rc1[p1] * Rcaa2[(q2*n_orb2 + s2)*n_orb2 + r2];
				}
			}
		}
	}
	H *= 2;
	for (PyInt p1=0;  p1<n_orb1;  p1++)
	{
		for (PyInt q2=0;  q2<n_orb2;  q2++)
		{
			H += h[p1*n_orb2 + q2] * Rc1[p1] * Ra2[q2];
		}
	}
	if(TIME_INSIDE_FUNCTION) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		printf("dimer_1min1pls finished execution in %f seconds\n", t);
	}
	return H;
}



Double dimer_00(PyInt n_orb1, PyInt n_orb2, Double* Rca1, Double* Rca2, Double* V)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	Double H = 0;
	for (PyInt p1=0;  p1<n_orb1;  p1++)
	{
		for (PyInt q2=0;  q2<n_orb2;  q2++)
		{
			for (PyInt r1=0;  r1<n_orb1;  r1++)
			{
				for (PyInt s2=0;  s2<n_orb2;  s2++)
				{
					int v_index = ((p1*n_orb2 + q2)*n_orb1 + r1)*n_orb2 + s2;
					int rca1_index = p1*n_orb1 +r1;
					int rca2_index = q2 * n_orb2 + s2;
					H += V[((p1*n_orb2 + q2)*n_orb1 + r1)*n_orb2 + s2] * Rca1[p1*n_orb1 + r1] * Rca2[q2*n_orb2 + s2];
				}
			}
		}
	}
	if(TIME_INSIDE_FUNCTION) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		printf("dimer_00 finished execution in %f seconds\n", t); 
	}
	return 4*H;
}



Double dimer_1pls1min(PyInt n_orb1, PyInt n_orb2, Double* Ra1, Double* Rcaa1, Double* Rc2, Double* Rcca2, Double* h, Double* V2221, Double* V2111)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	Double H = 0;
	for (PyInt p2=0;  p2<n_orb2;  p2++)
	{
		for (PyInt q2=0;  q2<n_orb2;  q2++)
		{
			for (PyInt r2=0;  r2<n_orb2;  r2++)
			{
				for (PyInt s1=0;  s1<n_orb1;  s1++)
				{
					H += V2221[((p2*n_orb2 + q2)*n_orb2 + r2)*n_orb1 + s1] * Rcca2[(q2*n_orb2 + p2)*n_orb2 + r2] * Ra1[s1];
				}
			}
		}
	}
	for (PyInt p2=0;  p2<n_orb2;  p2++)
	{
		for (PyInt q1=0;  q1<n_orb1;  q1++)
		{
			for (PyInt r1=0;  r1<n_orb1;  r1++)
			{
				for (PyInt s1=0;  s1<n_orb1;  s1++)
				{
					H += V2111[((p2*n_orb1 + q1)*n_orb1 + r1)*n_orb1 + s1] * Rc2[p2] * Rcaa1[(q1*n_orb1 + s1)*n_orb1 + r1];
				}
			}
		}
	}
	H *= 2;
	for (PyInt p2=0;  p2<n_orb2;  p2++)
	{
		for (PyInt q1=0;  q1<n_orb1;  q1++)
		{
			H += h[p2*n_orb1 + q1] * Rc2[p2] * Ra1[q1];
		}
	}
	if(TIME_INSIDE_FUNCTION) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		printf("dimer_1pls1min finished execution in %f seconds\n",t); 
	}
	return H;
	}



Double dimer_2pls2min(PyInt n_orb1, PyInt n_orb2, Double* Raa1, Double* Rcc2, Double* V)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	Double H = 0;
	for (PyInt p2=0;  p2<n_orb2;  p2++)
	{
		for (PyInt q2=0;  q2<n_orb2;  q2++)
		{
			for (PyInt r1=0;  r1<n_orb1;  r1++)
			{
				for (PyInt s1=0;  s1<n_orb1;  s1++)
				{
					H += V[((p2*n_orb2 + q2)*n_orb1 + r1)*n_orb1 + s1] * Rcc2[p2*n_orb2 + q2] * Raa1[s1*n_orb1 + r1];
				}
			}
		}
	}
	if(TIME_INSIDE_FUNCTION) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		printf("dimer_2pls2min finished execution in %f seconds\n", t);
	}
	return H;
}





void monomer_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb, Double** Rca,  Double** Rccaa, Double** h, Double** V)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	for (PyInt n=0;  n<n_elem;  n++)
	{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = monomer(n_orb[n], Rca[n], Rccaa[n], h[n], V[n]);
	}
	if(TIME_LOOP) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		total_time += t;
		printf("monomer_loop finished execution in %f seconds\n",t); 
	}
	return;
}

void dimer_2min2pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Rcc1, Double** Raa2, Double** V)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	for (PyInt n=0;  n<n_elem;  n++)
	{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = dimer_2min2pls(n_orb1[n], n_orb2[n], Rcc1[n], Raa2[n], V[n]);
	}
	if(TIME_LOOP) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		total_time += t;
		printf("dimer_2min2pls_loop finished execution in %f seconds\n", t);
	}
	return;
}

void dimer_1min1pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Rc1,  Double** Rcca1, Double** Ra2, Double** Rcaa2, Double** h, Double** V1, Double** V2)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	for (PyInt n=0;  n<n_elem;  n++)
	{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = sign[n] * dimer_1min1pls(n_orb1[n], n_orb2[n], Rc1[n], Rcca1[n], Ra2[n], Rcaa2[n], h[n], V1[n], V2[n]);
	}
	if(TIME_LOOP) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		total_time += t;
		printf("dimer_1min1pls_loop finished execution in %f seconds\n", t); 
	}
	return;
}

void dimer_00_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Rca1, Double** Rca2, Double** V)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	/*
	for(PyInt n=0; n<n_elem; n++)
	{
		PyInt index = i[n]*dim[n]+j[n];
		H[n][index] = dimer_00_cuda(n_orb1[n], n_orb2[n], Rca1[n], Rca2[n], V[n])
	}
	*/
	for (PyInt n=0;  n<n_elem;  n++)
	{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = dimer_00(n_orb1[n], n_orb2[n], Rca1[n], Rca2[n], V[n]);
	}
	if(TIME_LOOP) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		total_time += t;
		printf("dimer_00_loop finished execution in %f seconds\n", t); 
	}
	return;
}

void dimer_1pls1min_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Ra1,  Double** Rcaa1, Double** Rc2, Double** Rcca2, Double** h, Double** V1, Double** V2)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	for (PyInt n=0;  n<n_elem;  n++)
	{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = sign[n] * dimer_1pls1min(n_orb1[n], n_orb2[n], Ra1[n], Rcaa1[n], Rc2[n], Rcca2[n], h[n], V1[n], V2[n]);
	}
	if(TIME_LOOP) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		total_time += t;
		printf("dimer_1pls1min_loop finished execution in %f seconds\n", t); 
	}
	return;
}

void dimer_2pls2min_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Raa1, Double** Rcc2, Double** V)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	for (PyInt n=0;  n<n_elem;  n++)
	{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = dimer_2pls2min(n_orb1[n], n_orb2[n], Raa1[n], Rcc2[n], V[n]);
	}
	if(TIME_LOOP) {
		gettimeofday(&t2,NULL);
		double t = (double)(t2.tv_usec-t1.tv_usec)/1000000 + (double)(t2.tv_sec-t1.tv_sec);
		total_time += t;
		printf("dimer_2pls2min_loop finished execution in %f seconds\n", t);
	}
	return;
}

void print_total_time()
{
	printf("Total: %f seconds.\n", total_time);
}
