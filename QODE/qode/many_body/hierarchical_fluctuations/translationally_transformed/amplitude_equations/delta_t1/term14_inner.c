/*
 * header
 */

// This is automatically generated code.  Do not edit directly!

#include <stdint.h>

void loops(int64_t n_occ, int64_t n_vrt, double* O1, double* F, double* V, double* T1, double* T2, int64_t O1_0, int64_t F_0, int64_t V_0, int64_t T1_0, int64_t T2_0)
    {
    int64_t n_orb = n_occ + n_vrt;

    for (int64_t i=0; i<n_occ; i++)
        {
        int64_t O1_i = O1_0 + i;
        int64_t  F_i =  F_0;
        int64_t  V_i =  V_0;
        int64_t T1_i = T1_0 + i;
        int64_t T2_i = T2_0;
        for (int64_t a=0; a<n_vrt; a++)
            {
            int64_t O1_ai = O1_i + a*n_occ;
            int64_t  F_ai =  F_i;
            int64_t  V_ai =  V_i;
            int64_t T1_ai = T1_i;
            int64_t T2_ai = T2_i + a*n_occ*n_occ;
            for (int64_t k=0; k<n_occ; k++)
                {
                int64_t O1_aik = O1_ai;
                int64_t  F_aik =  F_ai;
                int64_t  V_aik =  V_ai + k*n_orb*n_orb*n_orb;
                int64_t T1_aik = T1_ai;
                int64_t T2_aik = T2_ai + k*n_occ;
                for (int64_t c=0; c<n_vrt; c++)
                    {
                    int64_t O1_aikc = O1_aik;
                    int64_t  F_aikc =  F_aik;
                    int64_t  V_aikc =  V_aik + (c+n_occ)*n_orb;
                    int64_t T1_aikc = T1_aik;
                    int64_t T2_aikc = T2_aik + c*n_vrt*n_occ*n_occ;
                    for (int64_t l=0; l<n_occ; l++)
                        {
                        int64_t O1_aikcl = O1_aikc;
                        int64_t  F_aikcl =  F_aikc;
                        int64_t  V_aikcl =  V_aikc + l*n_orb*n_orb;
                        int64_t T1_aikcl = T1_aikc;
                        int64_t T2_aikcl = T2_aikc + l;
                        for (int64_t d=0; d<n_vrt; d++)
                            {
                            int64_t O1_aikcld = O1_aikcl;
                            int64_t  F_aikcld =  F_aikcl;
                            int64_t  V_aikcld =  V_aikcl + (d+n_occ);
                            int64_t T1_aikcld = T1_aikcl + d*n_occ;
                            int64_t T2_aikcld = T2_aikcl;
                            O1[O1_aikcld] += -0.5 * V[V_aikcld] * T2[T2_aikcld] * T1[T1_aikcld];
                            }
                        }
                    }
                }
            }
        }

    return;
    }
