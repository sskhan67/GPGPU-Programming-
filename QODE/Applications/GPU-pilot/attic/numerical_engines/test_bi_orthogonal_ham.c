/*   (C) Copyright 2018 Yuhong Liu
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
/*

Two-body Hamiltonian contraction routine

For biorthogonal h and V matrices

< psi_I psi_J | \sum_{pq} h^p_q p+ q + \sum_{pqrs} v^{pq}_{rs} p+ q+ s r | psi_I' psi_J'>

= \sum_{pq} h^p_q < psi_I psi_J | p+ q | psi_I' psi_J'>  + 
  \sum_{pqrs} v^{pq}_{rs} < psi_I psi_J | p+ q+ s r | psi_I' psi_J'> 

h and V are semi-molecular orbitals for DIMER.

*/

#include <stdio.h>
#include <unistd.h> // for ssize_t : signed size_t 
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>


typedef struct 
{
	char   name[5];
	double *** data;
} Storage;

typedef struct
{
	ssize_t num_elec;
	ssize_t * pairs; // dimension = n_pairs x 2
	ssize_t n_pairs;
	double  ** data;
}
ContractedHamiltonian;


enum operators     // The order of operator sequences is fixed
{                  // The numbers here are acutally the indices for the Storage array
	seq_a    = 0,
	seq_c    = 1,
	seq_aa   = 2,
	seq_cc   = 3,
	seq_ca   = 4,
	seq_caa  = 5,
	seq_cca  = 6,
	seq_ccaa = 7
};

void show_storage(Storage store[])
{
	// Always an array of 8 elements
	for (int i = 0; i < 8; i++)
	{
		printf("Storage[%d] name = %4s;  data = %p\n", i, store[i].name, store[i].data);
	}
}

void show_cont_ham(ContractedHamiltonian contd_H[], ssize_t dim)
{
	for (ssize_t i = 0; i < dim; i++)
	{
		printf("H[%zd] --> num elec = %zd num pairs = %zd data ptr = %p\n", i, contd_H[i].num_elec, contd_H[i].n_pairs, contd_H[i].data);
		for (ssize_t x = 0; x < contd_H[i].n_pairs; x++)
		{
			printf("index_list[%zd] = [%zd,%zd]\n", x, contd_H[i].pairs[x*2], contd_H[i].pairs[x*2+1]);
		}
	}
}

ssize_t charge_states_to_contH_index(ssize_t ch_st_index1, ssize_t ch_st_index2, ssize_t num_elec1[], ssize_t num_elec2[], 
                                     ContractedHamiltonian contd_H[], ssize_t num_H_mat)
{
	ssize_t ind = 0;
	ssize_t num_elec = num_elec1[ch_st_index1] + num_elec2[ch_st_index2];
	while (ind < num_H_mat)
	{
		if (num_elec == contd_H[ind].num_elec)
		{
			return ind;
		}
		ind++;
	}
	puts("ERROR in function charge_states_to_contH_index");
	return -1;
}

ssize_t charge_states_to_contH_data_index(ssize_t ch_st_index1, ssize_t ch_st_index2, ssize_t ch_st_index3, ssize_t ch_st_index4, // Hamiltonian order <1,2|H|3,4>
	                                      ContractedHamiltonian contd_H[], ssize_t data_index)
{
	ssize_t * index_list = contd_H[data_index].pairs;
	ssize_t loc1 = -1, loc2 = -1;
	for (ssize_t i = 0; i < contd_H[data_index].n_pairs; i++)
	{
		if ( index_list[i*2] == ch_st_index1 && index_list[i*2+1] == ch_st_index2 ) { loc1 = i; }
		if ( index_list[i*2] == ch_st_index3 && index_list[i*2+1] == ch_st_index4 ) { loc2 = i; }
	}

	if (loc1 == -1 || loc2 == -1) { puts("ERROR in function charge_states_to_contH_data_index"); }

	assert(loc1 != -1);
	assert(loc2 != -1);
	return loc1 * contd_H[data_index].n_pairs + loc2;
}

typedef enum 
{
	frag1,
	frag2
} Fragment;

ssize_t monomer_spin_index_to_dimer_spatial(ssize_t spin_idx, ssize_t n_spat1, ssize_t n_spat2, Fragment frag)
{
	if (frag == frag1)
	{
		if (spin_idx < n_spat1) { return spin_idx; }
		else { return spin_idx - n_spat1; }
	}
	else
	{
		if (spin_idx < n_spat2) { return spin_idx + n_spat1; }
		else { return spin_idx - n_spat2 + n_spat1; }
	}
}

void check_bounds(ssize_t index, ssize_t bound)
{
	if (index < bound) {}
	else { printf("%zd/%zd out of bound!\n", index, bound); }
}



void contract_biorthogonal_ham(ssize_t n_charge_st1, ssize_t n_spin1, ssize_t num_states_used1[], ssize_t num_elec1[],
	                           ssize_t n_charge_st2, ssize_t n_spin2, ssize_t num_states_used2[], ssize_t num_elec2[],
	                           Storage storage1[], Storage storage2[],
	                           double h_dimer[], double V_dimer[],  // both in spatial orbitals 
	                           ContractedHamiltonian contd_H[], ssize_t num_H_mat, ssize_t nthd)
{
	ssize_t n_spat1 = n_spin1 / 2;   // num spatial orbs for frag1
	ssize_t n_spat2 = n_spin2 / 2;
	ssize_t n_spat_dimer = n_spat1 + n_spat2;  // num spatial orbs for two body
	ssize_t pp, qq, rr, ss;
	ssize_t H_index, h_index, V_index; 
	ssize_t ten1_index, ten2_index;
	ssize_t index1, index2, mat_index;
	double *** p_data1   = NULL, *** p_data2   = NULL;
	double *   p_tensor1 = NULL, *   p_tensor2 = NULL;
	double * H_mat;

	// omp_set_dynamic(0);
    // omp_set_num_threads(nthd);

	// Eqn: <I|<K| H |J>|L> --> \sum ... <I|b1|J> <K|b2|L> 

	// #pragma omp for
	for (ssize_t I = 0; I < n_charge_st1; I++)  
	{
		for (ssize_t J = 0; J < n_charge_st1; J++)
		{
			for (ssize_t K = 0; K < n_charge_st2; K++)
			{
				for (ssize_t L = 0; L < n_charge_st2; L++)
				{
					index1 = charge_states_to_contH_index(I, K, num_elec1, num_elec2, contd_H, num_H_mat);
					index2 = charge_states_to_contH_index(J, L, num_elec1, num_elec2, contd_H, num_H_mat);
					if (index1 == index2)
					{
					

						// Two operators: (4 cases)

						/*
						// 1)  I_{A}  c_{B}(p) a_{B}(q) = c_{B}(p) a_{B}(q)

						p_data1 = NULL; 
						p_data2 = storage2[seq_ca].data;
						p_tensor1 = NULL;

						if (I == J)  // frag1 same state
						{
							if ( p_data2[K*n_charge_st2+L] ) 
							{
								mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
								H_mat     = contd_H[index1].data[mat_index];
								for (ssize_t i = 0; i < num_states_used1[I]; i++)
								{
									for (ssize_t j = 0; j < num_states_used1[J]; j++)
									{
										if (i == j)
										{
											for (ssize_t k = 0; k < num_states_used2[K]; k++)
											{
												for (ssize_t l = 0; l < num_states_used2[L]; l++)
												{

													p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
													for (size_t p = 0; p < n_spin2; p++)
													{
														pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag2);
														for (size_t q = 0; q < n_spin2; q++)
														{
															qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag2);
															H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
															h_index = pp * n_spat_dimer + qq;
															ten2_index = p*n_spin2+q;
															// #pragma omp atomic
															// H_mat[ H_index ] += h_dimer[h_index] * p_tensor2[ten2_index];  
														}
													}
												}
											}
										}
									}
								}

							}
						}
						*/


						// 2)  a_{A}(q) c_{B}(p) = - c_{B}(p) a_{A}(q)

						p_data1 = storage1[seq_a].data; 
						p_data2 = storage2[seq_c].data;
						if ( p_data1[I*n_charge_st1+J] && p_data2[K*n_charge_st2+L] ) 
						{
							mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
							H_mat     = contd_H[index1].data[mat_index];
							for (ssize_t i = 0; i < num_states_used1[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used1[J]; j++)
								{
									p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
									for (ssize_t k = 0; k < num_states_used2[K]; k++)
									{
										for (ssize_t l = 0; l < num_states_used2[L]; l++)
										{

											p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
											for (size_t p = 0; p < n_spin2; p++)
											{
												pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag2);
												for (size_t q = 0; q < n_spin1; q++)
												{
													qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag1);
													H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
													h_index = pp * n_spat_dimer + qq;
													// #pragma omp atomic
													printf("h value = %.3e\n", -h_dimer[h_index]);
													H_mat[H_index] += -h_dimer[h_index] * p_tensor1[q] * p_tensor2[p];
												}
											}
										}
									}

								}
							}
						}



						// 3) c_{A}(p) a_{B}(q) = c_{A}(p) a_{B}(q)

						p_data1 = storage1[seq_c].data; 
						p_data2 = storage2[seq_a].data;
						if ( p_data1[I*n_charge_st1+J] && p_data2[K*n_charge_st2+L] ) 
						{
							mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
							H_mat     = contd_H[index1].data[mat_index];
							for (ssize_t i = 0; i < num_states_used1[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used1[J]; j++)
								{
									p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
									for (ssize_t k = 0; k < num_states_used2[K]; k++)
									{
										for (ssize_t l = 0; l < num_states_used2[L]; l++)
										{

											p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
											for (size_t p = 0; p < n_spin1; p++)
											{
												pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag1);
												for (size_t q = 0; q < n_spin2; q++)
												{
													qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag2);
													H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
													h_index = pp * n_spat_dimer + qq;
													// #pragma omp atomic
													printf("h value = %.3e\n", h_dimer[h_index]);
													H_mat[H_index] += h_dimer[h_index] * p_tensor1[p] * p_tensor2[q];
												}
											}
										}
									}

								}
							}
						}

						/*
						// 4)  c_{A}(p) a_{A}(q)   I_{B} = c_{A}(p) a_{A}(q)

						p_data1 = storage1[seq_ca].data; 
						p_data2 = NULL;

						if (K == L)  // frag1 same state 
						{
							if ( p_data1[I*n_charge_st1+J] )
							{
								mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
								H_mat     = contd_H[index1].data[mat_index];
								for (ssize_t i = 0; i < num_states_used1[I]; i++)
								{
									for (ssize_t j = 0; j < num_states_used1[J]; j++)
									{
										if (i == j)
										{
											p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
											for (ssize_t k = 0; k < num_states_used2[K]; k++)
											{
												for (ssize_t l = 0; l < num_states_used2[L]; l++)
												{
													for (size_t p = 0; p < n_spin1; p++)
													{
														pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag1);
														for (size_t q = 0; q < n_spin1; q++)
														{
															qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag1);
															H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
															h_index = pp * n_spat_dimer + qq;
															ten1_index = p * n_spin1 + q;
															// #pragma omp atomic
															// H_mat[H_index] += h_dimer[h_index] * p_tensor1[ten1_index];
														}
													}
												}
											}
										}
									}
								}

							}
						}
						*/
						
						// Four operators: (9 cases)

						/*
						// 1) I_{A}  c_{B1}(p) c_{B2}(q) a_{B1}(s) a_{B2}(r) = c_{B1}(p) c_{B2}(q) a_{B1}(s) a_{B2}(r) 

						p_data1 = NULL;
						p_data2 = storage2[seq_ccaa].data;
						p_tensor1 = NULL;

						
						if (I == J)  // frag1 same state
						{
							if ( p_data2[K*n_charge_st2+L] ) 
							{
								mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
								H_mat     = contd_H[index1].data[mat_index];
								for (ssize_t i = 0; i < num_states_used1[I]; i++)
								{
									for (ssize_t j = 0; j < num_states_used1[J]; j++)
									{
										if (i == j)
										{
											for (ssize_t k = 0; k < num_states_used2[K]; k++)
											{
												for (ssize_t l = 0; l < num_states_used2[L]; l++)
												{
													p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
													for (size_t p = 0; p < n_spin2; p++)
													{
														pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag2);
														for (size_t q = 0; q < n_spin2; q++)
														{
															qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag2);
															for (size_t r = 0; r < n_spin2; r++)
															{
																rr = monomer_spin_index_to_dimer_spatial(r, n_spat1, n_spat2, frag2);
																for (size_t s = 0; s < n_spin2; s++)
																{
																	ss = monomer_spin_index_to_dimer_spatial(s, n_spat1, n_spat2, frag2);
																	H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
																	V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
																	ten2_index = (p * n_spin2 + q) * n_spin2 * n_spin2 + s * n_spin2 + r;
																	// #pragma omp atomic
																	// H_mat[H_index] += V_dimer[V_index] * p_tensor2[ten2_index];  
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
						*/

						// 2)  a_{A}(s)   c_{B1}(p) c_{B2}(q) a_{B}(r) =      c_{B1}(p) c_{B2}(q) a_{A}(s) a_{B}(r)
						//                                             = -1 * c_{B1}(p) c_{B2}(q) a_{B}(r) a_{A}(s)
						p_data1 = storage1[seq_a].data; 
						p_data2 = storage2[seq_cca].data;
						if ( p_data1[I*n_charge_st1+J] && p_data2[K*n_charge_st2+L] ) 
						{
							mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
							H_mat     = contd_H[index1].data[mat_index];
							for (ssize_t i = 0; i < num_states_used1[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used1[J]; j++)
								{
									p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
									for (ssize_t k = 0; k < num_states_used2[K]; k++)
									{
										for (ssize_t l = 0; l < num_states_used2[L]; l++)
										{
											p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
											for (size_t p = 0; p < n_spin2; p++)
											{
												pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag2);
												for (size_t q = 0; q < n_spin2; q++)
												{
													qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag2);
													for (size_t r = 0; r < n_spin2; r++)
													{
														rr = monomer_spin_index_to_dimer_spatial(r, n_spat1, n_spat2, frag2);
														for (size_t s = 0; s < n_spin1; s++)
														{
															ss = monomer_spin_index_to_dimer_spatial(s, n_spat1, n_spat2, frag1);
															ten2_index = (p * n_spin2 + q) * n_spin2  + r;
															H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
															
															// first eqn
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															// #pragma omp atomic
															printf("V value = %.3e\n", V_dimer[V_index]);
															H_mat[H_index] += V_dimer[V_index] * p_tensor1[s] * p_tensor2[ten2_index];  
															// second eqn
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + ss * n_spat_dimer + rr;
															// #pragma omp atomic
															printf("V value = %.3e\n", -V_dimer[V_index]);
															H_mat[H_index] += -V_dimer[V_index] * p_tensor1[s] * p_tensor2[ten2_index]; 
														}
													}
												}
											}
										}
									}
								}
							}		
						}

						// 3) c_{A}(p) c_{B}(q) a_{B1}(s) a_{B2}(r) =       c_{A}(p) c_{B}(q) a_{B1}(s) a_{B2}(r)
						//                                          =  -1 * c_{B}(q) c_{A}(p) a_{B1}(s) a_{B2}(r)

						p_data1 = storage1[seq_c].data; 
						p_data2 = storage2[seq_caa].data;
						if ( p_data1[I*n_charge_st1+J] && p_data2[K*n_charge_st2+L] ) 
						{
							mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
							H_mat     = contd_H[index1].data[mat_index];
							for (ssize_t i = 0; i < num_states_used1[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used1[J]; j++)
								{
									p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
									for (ssize_t k = 0; k < num_states_used2[K]; k++)
									{
										for (ssize_t l = 0; l < num_states_used2[L]; l++)
										{
											p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
											for (size_t p = 0; p < n_spin1; p++)
											{
												pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag1);
												for (size_t q = 0; q < n_spin2; q++)
												{
													qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag2);
													for (size_t r = 0; r < n_spin2; r++)
													{
														rr = monomer_spin_index_to_dimer_spatial(r, n_spat1, n_spat2, frag2);
														for (size_t s = 0; s < n_spin2; s++)
														{
															ss = monomer_spin_index_to_dimer_spatial(s, n_spat1, n_spat2, frag2);
															ten2_index = (q * n_spin2 + s) * n_spin2 + r;
															H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
															// first eqn
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															// #pragma omp atomic
															printf("V value = %.3e\n", V_dimer[V_index]);
															H_mat[H_index] += V_dimer[V_index] * p_tensor1[p] * p_tensor2[ten2_index];  
															
															// second eqn
															V_index = (qq * n_spat_dimer + pp) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															// #pragma omp atomic
															printf("V value = %.3e\n", -V_dimer[V_index]);
															H_mat[H_index] += -V_dimer[V_index] * p_tensor1[p] * p_tensor2[ten2_index];  
														}
													}
												}
											}
										}
									}
								}
							}
						}



						// 4) a_{A1}(s) a_{A2}(r) c_{B1}(p) c_{B2}(q) = c_{B1}(p) c_{B2}(q) a_{A1}(s) a_{A2}(r) 

						p_data1 = storage1[seq_aa].data; 
						p_data2 = storage2[seq_cc].data;
						if ( p_data1[I*n_charge_st1+J] && p_data2[K*n_charge_st2+L] ) 
						{
							mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
							H_mat     = contd_H[index1].data[mat_index];
							for (ssize_t i = 0; i < num_states_used1[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used1[J]; j++)
								{
									p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
									for (ssize_t k = 0; k < num_states_used2[K]; k++)
									{
										for (ssize_t l = 0; l < num_states_used2[L]; l++)
										{
											p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
											for (size_t p = 0; p < n_spin2; p++)
											{
												pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag2);
												for (size_t q = 0; q < n_spin2; q++)
												{
													qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag2);
													for (size_t r = 0; r < n_spin1; r++)
													{
														rr = monomer_spin_index_to_dimer_spatial(r, n_spat1, n_spat2, frag1);
														for (size_t s = 0; s < n_spin1; s++)
														{
															ss = monomer_spin_index_to_dimer_spatial(s, n_spat1, n_spat2, frag1);
															H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															ten1_index = s * n_spin1 + r;
															ten2_index = p * n_spin2 + q;
															// #pragma omp atomic
															printf("V value = %.3e\n", V_dimer[V_index]);
															H_mat[H_index] += V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[ten2_index];  
														}
													}
												}
											}
										}
									}
								}
							}
						}



						// 5)  c_{A1}(p) c_{A2}(q) a_{B1}(s) a_{B2}(r) = c_{A1}(p) c_{A2}(q) a_{B1}(s) a_{B2}(r)
						p_data1 = storage1[seq_cc].data; 
						p_data2 = storage2[seq_aa].data;
						if ( p_data1[I*n_charge_st1+J] && p_data2[K*n_charge_st2+L] ) 
						{
							mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
							H_mat     = contd_H[index1].data[mat_index];
							for (ssize_t i = 0; i < num_states_used1[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used1[J]; j++)
								{
									p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
									for (ssize_t k = 0; k < num_states_used2[K]; k++)
									{
										for (ssize_t l = 0; l < num_states_used2[L]; l++)
										{
											p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
											for (size_t p = 0; p < n_spin1; p++)
											{
												pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag1);
												for (size_t q = 0; q < n_spin1; q++)
												{
													qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag1);
													for (size_t r = 0; r < n_spin2; r++)
													{
														rr = monomer_spin_index_to_dimer_spatial(r, n_spat1, n_spat2, frag2);
														for (size_t s = 0; s < n_spin2; s++)
														{
															ss = monomer_spin_index_to_dimer_spatial(s, n_spat1, n_spat2, frag2);
															H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															ten1_index = p * n_spin1 + q;
															ten2_index = s * n_spin2 + r;
															// #pragma omp atomic
															printf("V value = %.3e\n", V_dimer[V_index]);
															H_mat[H_index] += V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[ten2_index];  
														}
													}
												}
											}
										}
									}
								}
							}
						}



						// 6)   c_{A}(p) a_{A}(r) c_{B}(q) a_{B}(s) =      c_{A}(p) c_{B}(q) a_{B}(s) a_{A}(r)
						//                                          = -1 * c_{B}(q) c_{A}(p) a_{B}(s) a_{A}(r)
						//                                          = -1 * c_{A}(p) c_{B}(q) a_{A}(r) a_{B}(s)
						//                                          =      c_{B}(q) c_{A}(p) a_{A}(r) a_{B}(s)

						p_data1 = storage1[seq_ca].data; 
						p_data2 = storage2[seq_ca].data;
						if ( p_data1[I*n_charge_st1+J] && p_data2[K*n_charge_st2+L] ) 
						{
							mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
							H_mat     = contd_H[index1].data[mat_index];
							for (ssize_t i = 0; i < num_states_used1[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used1[J]; j++)
								{
									p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
									for (ssize_t k = 0; k < num_states_used2[K]; k++)
									{
										for (ssize_t l = 0; l < num_states_used2[L]; l++)
										{
											p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
											for (size_t p = 0; p < n_spin1; p++)
											{
												pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag1);
												for (size_t q = 0; q < n_spin2; q++)
												{
													qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag2);
													for (size_t r = 0; r < n_spin1; r++)
													{
														rr = monomer_spin_index_to_dimer_spatial(r, n_spat1, n_spat2, frag1);
														for (size_t s = 0; s < n_spin2; s++)
														{
															ss = monomer_spin_index_to_dimer_spatial(s, n_spat1, n_spat2, frag2);
															H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
															ten1_index = p * n_spin1 + r;
															ten2_index = q * n_spin2 + s;
															// 1st eqn
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															// #pragma omp atomic
															printf("V value = %.3e\n", V_dimer[V_index]);
															H_mat[H_index] += V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[ten2_index];  
															
															// 2nc eqn
															V_index = (qq * n_spat_dimer + pp) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															// #pragma omp atomic
															printf("V value = %.3e\n", -V_dimer[V_index]);
															H_mat[H_index] += -V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[ten2_index];  
															
															// 3rd eqn
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + ss * n_spat_dimer + rr;
															// #pragma omp atomic
															printf("V value = %.3e\n", -V_dimer[V_index]);
															H_mat[H_index] += -V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[ten2_index];  
															
															// 4th eqn
															V_index = (qq * n_spat_dimer + pp) * n_spat_dimer * n_spat_dimer + ss * n_spat_dimer + rr;
															// #pragma omp atomic
															printf("V value = %.3e\n", V_dimer[V_index]);
															H_mat[H_index] += V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[ten2_index];  
														}
													}
												}
											}
										}
									}
								}
							}
						}

						// 7) c_{A}(p) a_{A1}(s) a_{A2}(r) c_{B}(q) =      c_{A}(p) c_{B}(q) a_{A1}(s) a_{A2}(r) 
						//                                          = -1 * c_{B}(q) c_{A}(p) a_{A1}(s) a_{A2}(r)

						p_data1 = storage1[seq_caa].data; 
						p_data2 = storage2[seq_c].data;
						if ( p_data1[I*n_charge_st1+J] && p_data2[K*n_charge_st2+L] ) 
						{
							mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
							H_mat     = contd_H[index1].data[mat_index];
							for (ssize_t i = 0; i < num_states_used1[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used1[J]; j++)
								{
									p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
									for (ssize_t k = 0; k < num_states_used2[K]; k++)
									{
										for (ssize_t l = 0; l < num_states_used2[L]; l++)
										{
											p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
											for (size_t p = 0; p < n_spin1; p++)
											{
												pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag1);
												for (size_t q = 0; q < n_spin2; q++)
												{
													qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag2);
													for (size_t r = 0; r < n_spin1; r++)
													{
														rr = monomer_spin_index_to_dimer_spatial(r, n_spat1, n_spat2, frag1);
														for (size_t s = 0; s < n_spin1; s++)
														{
															ss = monomer_spin_index_to_dimer_spatial(s, n_spat1, n_spat2, frag1);
															ten1_index = (p * n_spin1 + s) * n_spin1 + r;
															H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
															// 1st eqn
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															// #pragma omp atomic
															printf("V value = %.3e\n", V_dimer[V_index]);
															H_mat[H_index] += V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[q];
															
															// 2nd eqn
															V_index = (qq * n_spat_dimer + pp) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															// #pragma omp atomic
															printf("V value = %.3e\n", -V_dimer[V_index]);
															H_mat[H_index] += -V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[q];
														}
													}
												}
											}
										}
									}
								}
							}
						}

						// 8) c_{A1}(p) c_{A2}(q) a_{A}(s) a_{B}(r) =      c_{A1}(p) c_{A2}(q) a_{A}(s) a_{B}(r) 
						//                                          = -1 * c_{A1}(p) c_{A2}(q) a_{B}(r) a_{A}(s) 

						p_data1 = storage1[seq_cca].data; 
						p_data2 = storage2[seq_a].data;
						if ( p_data1[I*n_charge_st1+J] && p_data2[K*n_charge_st2+L] ) 
						{
							mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
							H_mat     = contd_H[index1].data[mat_index];
							for (ssize_t i = 0; i < num_states_used1[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used1[J]; j++)
								{
									p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
									for (ssize_t k = 0; k < num_states_used2[K]; k++)
									{
										for (ssize_t l = 0; l < num_states_used2[L]; l++)
										{
											p_tensor2 = p_data2[K*n_charge_st2+L][k*num_states_used2[L]+l];
											for (size_t p = 0; p < n_spin1; p++)
											{
												pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag1);
												for (size_t q = 0; q < n_spin1; q++)
												{
													qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag1);
													for (size_t r = 0; r < n_spin2; r++)
													{
														rr = monomer_spin_index_to_dimer_spatial(r, n_spat1, n_spat2, frag2);
														for (size_t s = 0; s < n_spin1; s++)
														{
															ss = monomer_spin_index_to_dimer_spatial(s, n_spat1, n_spat2, frag1);
															ten1_index = (p * n_spin1 + q) * n_spin1  + s;
															H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
															// 1st eqn
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
															// #pragma omp atomic
															printf("V value = %.3e\n", V_dimer[V_index]);
															H_mat[H_index] += V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[r];  
															
															// 2nd eqn
															V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + ss * n_spat_dimer + rr;
															// #pragma omp atomic
															printf("V value = %.3e\n", -V_dimer[V_index]);
															H_mat[H_index] += -V_dimer[V_index] * p_tensor1[ten1_index] * p_tensor2[r]; 
														}
													}
												}
											}
										}
									}
								}
							}
						}

						/*
						// 9)  c_{A1}(p) c_{A2}(q) a_{A1}(s) a_{A2}(r) = c_{A1}(p) c_{A2}(q) a_{A1}(s) a_{A2}(r)
						
						p_data1 = storage1[seq_ccaa].data; 
						p_data2 = NULL;

						if (K == L)  // frag1 same state
						{
							if ( p_data1[I*n_charge_st1+J] ) 
							{
								mat_index = charge_states_to_contH_data_index(I, K, J, L, contd_H, index1);
								H_mat     = contd_H[index1].data[mat_index];
								for (ssize_t i = 0; i < num_states_used1[I]; i++)
								{
									for (ssize_t j = 0; j < num_states_used1[J]; j++)
									{
										p_tensor1 = p_data1[I*n_charge_st1+J][i*num_states_used1[J]+j];
										for (ssize_t k = 0; k < num_states_used2[K]; k++)
										{
											for (ssize_t l = 0; l < num_states_used2[L]; l++)
											{
												if (k == l)
												{
													for (size_t p = 0; p < n_spin1; p++)
													{
														pp = monomer_spin_index_to_dimer_spatial(p, n_spat1, n_spat2, frag1);
														for (size_t q = 0; q < n_spin1; q++)
														{
															qq = monomer_spin_index_to_dimer_spatial(q, n_spat1, n_spat2, frag1);
															for (size_t r = 0; r < n_spin1; r++)
															{
																rr = monomer_spin_index_to_dimer_spatial(r, n_spat1, n_spat2, frag1);
																for (size_t s = 0; s < n_spin1; s++)
																{
																	ss = monomer_spin_index_to_dimer_spatial(s, n_spat1, n_spat2, frag1);
																	H_index = (i * num_states_used2[K] + k) * num_states_used1[J] * num_states_used2[L] + j * num_states_used2[L] + l;
																	V_index = (pp * n_spat_dimer + qq) * n_spat_dimer * n_spat_dimer + rr * n_spat_dimer + ss;
																	ten1_index = (p * n_spin1 + q) * n_spin1 * n_spin1 + s * n_spin1 + r;
																	// check_bounds(V_index, n_spat_dimer*n_spat_dimer*n_spat_dimer*n_spat_dimer);
																	// #pragma omp atomic
																	// H_mat[H_index] += V_dimer[V_index] * p_tensor1[ten1_index];  
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
						*/


					}  // index1 == index2
				}  // loop L
			}  // loop K
		}  //  loop J
	}  //  loop I




} // int main {}






















