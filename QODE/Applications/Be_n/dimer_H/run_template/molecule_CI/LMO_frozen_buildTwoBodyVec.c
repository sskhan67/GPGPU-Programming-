/*   (C) Copyright 2018 Anthony D. Dutoi
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
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <math.h>
#include "omp.h"
int compare (const void* a, const void* b) {return ( *(long long*)a - *(long long*)b );}

//
//  This turns A and B vectors into dimer FCI vector
//



double get_sign(long long* unsorted_config, long long dim)
{
	int sign = 1;
	for(long long i=0; i<dim; i++)
	{
		for(long long j=0; j<i; j++)
		{
			if (unsorted_config[j] > unsorted_config[i])
			{
				sign *= -1;
			}
		}
	}
	return (double)sign;
}

long long configToFciIndex(long long* config, long long num_elec, long long num_spin_orbs, long long* combo_matrix)
{
    long long index_num = 0;
    for(long long n=1; n<config[0]+1; n++)
    {
        index_num += combo_matrix[num_elec*num_spin_orbs - n];
    }
    for(long long i=0; i<num_elec-1; i++)
    {
        for(long long n=config[i]+2; n<config[i+1]+1; n++)
        {
            index_num += combo_matrix[(num_elec-i-1)*num_spin_orbs - n];
        }
    }
    return index_num;
}

int no_pauli_exclusion(long long* configA, long long* configB, long long num_elec_A, long long num_elec_B)
{
	for(long long i=0; i<num_elec_A; i++)
	{
		for(long long j=0; j<num_elec_B; j++)
		{
			if ( configA[i] == configB[j] )
			{
				return 0;	// False:  this DOES violate the Pauli exclusion principle
			}
		}
	}
	return 1;
}



void Transform_ABvecs_Into_TwoBody_Basis( double*    TwoBodyVec,
                                          long long  num_spin_orbs,
                                          long long  num_core_elec,
                                          long long  num_valence_elec_A,
                                          long long  num_valence_elec_B,
                                          double*    atomA_state,
                                          double*    atomB_state,
                                          long long* atomA_configs,
                                          long long* atomB_configs,
                                          long long  num_atomA_configs,
                                          long long  num_atomB_configs,
                                          long long* comboMat )
{
	long long num_valence_elec  = num_valence_elec_A + num_valence_elec_B;
	long long num_valence_orbs  = num_spin_orbs - num_core_elec;
	long long num_spat_orbs     = num_spin_orbs / 2;
	long long num_core_spat_orb = num_core_elec / 2;

	/*
	printf("TwoBodyVec         = %p \n",TwoBodyVec);
	printf("num_spin_orbs      = %d \n",num_spin_orbs);
	printf("num_core_elec      = %d \n",num_core_elec);
	printf("num_valence_elec_A = %d \n",num_valence_elec_A);
	printf("num_valence_elec_B = %d \n",num_valence_elec_B);
	printf("atomA_state        = %p \n",atomA_state);
	printf("atomB_state        = %p \n",atomB_state);
	printf("atomA_configs      = %p \n",atomA_configs);
	printf("atomB_configs      = %p \n",atomB_configs);
	printf("num_atomA_configs  = %d \n",num_atomA_configs);
	printf("num_atomB_configs  = %d \n",num_atomB_configs);
	printf("comboMat           = %p \n",comboMat);
	printf("num_valence_elec   = %d \n",num_valence_elec);
	printf("num_valence_orbs   = %d \n",num_valence_orbs);
	printf("num_spat_orbs      = %d \n",num_spat_orbs);
	printf("num_core_spat_orb  = %d \n",num_core_spat_orb);
	*/

	#pragma omp parallel for
	for(long long i=0; i<num_atomA_configs; i++)
	{
		for(long long j=0; j<num_atomB_configs; j++)
		{
			if ( no_pauli_exclusion( &atomA_configs[i*num_valence_elec_A ], &atomB_configs[j*num_valence_elec_B], \
				                                  num_valence_elec_A,                    num_valence_elec_B ) )
			{
				long long dimer_config[num_valence_elec];	// must be inside loop because of OMP

				for (long long k=0; k<num_valence_elec_A; k++) {dimer_config[k]                    = atomA_configs[i*num_valence_elec_A+k];}
				for (long long k=0; k<num_valence_elec_B; k++) {dimer_config[k+num_valence_elec_A] = atomB_configs[j*num_valence_elec_B+k];}

				for (long long k=0; k<num_valence_elec; k++)
				{
					if (dimer_config[k] > num_spat_orbs) {dimer_config[k] -= 2*num_core_spat_orb;}
					else                                 {dimer_config[k] -=   num_core_spat_orb;}
				}

				double sign = get_sign(dimer_config, num_valence_elec);
				qsort( dimer_config, num_valence_elec, sizeof(long long), compare );

				long long FCIindex = configToFciIndex( dimer_config, num_valence_elec, num_valence_orbs, comboMat );

				#pragma omp atomic
				TwoBodyVec[ FCIindex ] += sign * atomA_state[i] * atomB_state[j];
			} // pauli allowed
		} // Atom B
	} // Atom A

	// printf("did not make it here.\n");
	return;
}
