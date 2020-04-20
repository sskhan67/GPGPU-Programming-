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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> // for ssize_t : signed size_t



/* helper functions */

ssize_t value_in_array(ssize_t * array, ssize_t dim, ssize_t value)
{
	for (ssize_t i = 0; i < dim; i++)
	{ if (array[i] == value) {return 1;} }
	return 0;
}

ssize_t find_max(ssize_t array[], ssize_t dim)
{
	ssize_t max_val = array[0];
	for (ssize_t i = 1; i < dim; i++)
	{
		if ( max_val < array[i] )
		{
			max_val = array[i];
		}
		// else do not do anything
	}
	return max_val;
}

ssize_t find_min(ssize_t array[], ssize_t dim)
{
	ssize_t min_val = array[0];
	for (ssize_t i = 1; i < dim; i++)
	{
		if ( min_val > array[i] )
		{
			min_val = array[i];
		}
		// else do not do anything
	}
	return min_val;
}

typedef enum ValenceType {Core, Valence} ValType;

ValType get_val_type(ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t op_idx)
{
	if ( (op_idx >= 0 && op_idx <num_core_elec_alpha) || (op_idx >= num_spat_orbs && op_idx < num_spat_orbs + num_core_elec_alpha) )
	// alpha core || beta core
	{ return Core;}
	else
	{ return Valence; }
}

ssize_t find_index(ssize_t * array, ssize_t dim, ssize_t value)
{
	for (ssize_t i = 0; i < dim; i++)
	{
		if (array[i] == value) {return i;}
	}
	return -1; // not found
}

ssize_t total_index_to_valence_only_index(ssize_t index, ssize_t num_spat_orbs, ssize_t num_core_elec_alpha)
{  // num_core_elec_alpha == num_core_elec_beta
    if (index < num_spat_orbs) { return index - num_core_elec_alpha;     }
    else                       { return index - num_core_elec_alpha * 2; }
}

ssize_t find_FCI_valence_config_index(ssize_t num_val_elec, ssize_t num_val_spin_orbs, ssize_t* sorted_val_state, ssize_t* combo_mat)
{
	// sorted valence states is not a direct state from state1 or state2. All orbitals are decreased by
	// subtracting number of valence electrons to allow the configuration starts with 0,1,2,...
    ssize_t index_num = 0;
    for(ssize_t n=1; n<sorted_val_state[0]+1; n++)
    {
        index_num += combo_mat[num_val_elec * num_val_spin_orbs - n];
    }
    for(ssize_t i=0; i<num_val_elec-1; i++)
    {
        for(ssize_t n=sorted_val_state[i]+2; n<sorted_val_state[i+1]+1; n++)
        {
            index_num += combo_mat[(num_val_elec - i - 1) * num_val_spin_orbs - n];
        }
    }
    return index_num;
}

ssize_t find_config_index(ssize_t * state, ssize_t dim, ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t* combo_mat)
{
	// dim == num_elec in this case
	ssize_t num_val_elec = dim - num_core_elec_alpha*2;
	ssize_t num_val_orbs = (num_spat_orbs - num_core_elec_alpha)*2;
	ssize_t val_state[num_val_elec];
	ssize_t val_idx = 0;
	for (ssize_t i = 0; i< dim; i++)
	{
		if ( get_val_type(num_spat_orbs, num_core_elec_alpha, state[i]) == Valence )
		{
			val_state[val_idx] = total_index_to_valence_only_index(state[i], num_spat_orbs, num_core_elec_alpha);
			val_idx++;
		}
	}
	return find_FCI_valence_config_index( num_val_elec, num_val_orbs, val_state, combo_mat );
}

int compare(const void * a, const void * b)
{
	return ( *(ssize_t*)a - *(ssize_t*)b );
}

ssize_t n_swap_get_sign(ssize_t n)
{
	if (n % 2 ==0) { return  1;}
	else           { return -1;}
}

void show_state(ssize_t state[], ssize_t dim)
{

	for (ssize_t i = 0; i < dim; i++)
	{printf("%2zd ", state[i]);}

}

/* end of helper functions */



/* operator "b" functions */

ssize_t smart_a(ssize_t * P, ssize_t p, ssize_t * state_Q, ssize_t dim,
                ssize_t min_num_elec, ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t * combo_mat)
{
	ssize_t n_creation = -1;
	ssize_t new_dim = dim + n_creation;
	if (new_dim < min_num_elec) {return 0;}            // too few electrons
	if (!value_in_array(state_Q, dim, p)) {return 0;}  // annihilate on virtual
	if (get_val_type(num_spat_orbs, num_core_elec_alpha, p) == Core) {return 0;} // annihilation on frozen core is zero
	ssize_t new_state_Q[new_dim];
	ssize_t idx = 0;
	for (ssize_t i = 0; i < dim; i++) { if (state_Q[i] != p) { new_state_Q[idx++] = state_Q[i]; } }
	*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat); // put new config index in P
	ssize_t p_pos = find_index(state_Q, dim, p);
	return n_swap_get_sign(p_pos); // n swaps to move orb p to the left side of state Q
}

ssize_t smart_c(ssize_t * P,  ssize_t p, ssize_t * state_Q, ssize_t dim,
                ssize_t max_num_elec, ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t * combo_mat)
{
	ssize_t n_creation = 1;
	ssize_t new_dim = dim + n_creation;
	if (new_dim > max_num_elec) {return 0;}            // too many electrons
	if (value_in_array(state_Q, dim, p)) {return 0;}   // creation on occupied
	ssize_t new_state_Q[new_dim];
	for (ssize_t i = 0; i < dim; i++) {new_state_Q[i] = state_Q[i];}
	new_state_Q[new_dim - 1] = p;
	qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
	*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);  // put new config index in P
	ssize_t p_pos = find_index(new_state_Q, new_dim, p);
	return n_swap_get_sign(p_pos);
}

ssize_t smart_aa(ssize_t * P,  ssize_t p, ssize_t q, ssize_t * state_Q, ssize_t dim,
                 ssize_t min_num_elec, ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t * combo_mat)
{
	ssize_t n_creation = -2;
	ssize_t new_dim = dim + n_creation;
	if (new_dim < min_num_elec) {return 0;}   // too few electrons
	if (p == q) {return 0;}                   // annihilate same orbital two times
	if (!value_in_array(state_Q, dim, p)) {return 0;}  // annihilate a virtual
	if (!value_in_array(state_Q, dim, q)) {return 0;}  // annihilate a virtual
	if (get_val_type(num_spat_orbs, num_core_elec_alpha, p) == Core) {return 0;} // annihilation on frozen core is zero
	if (get_val_type(num_spat_orbs, num_core_elec_alpha, q) == Core) {return 0;} // annihilation on frozen core is zero
	ssize_t new_state_Q[new_dim];
	ssize_t idx = 0;
	for (ssize_t i = 0; i < dim; i++) { if ((state_Q[i] != p) && (state_Q[i] != q)) { new_state_Q[idx++] = state_Q[i]; } }
	*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
	ssize_t p_pos = find_index(state_Q, dim, p);
	ssize_t q_pos = find_index(state_Q, dim, q);
	if      (p_pos < q_pos) {return n_swap_get_sign(p_pos     + q_pos);}
	else if (p_pos > q_pos) {return n_swap_get_sign(p_pos - 1 + q_pos);}
	else                    {return 0;}
}

ssize_t smart_cc(ssize_t * P,  ssize_t p, ssize_t q, ssize_t * state_Q, ssize_t dim,
                 ssize_t max_num_elec, ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t * combo_mat)
{
	ssize_t n_creation = 2;
	ssize_t new_dim = dim + n_creation;
	if (new_dim > max_num_elec) {return 0;}  // too many electrons
	if (p == q) {return 0;}                  // create same orbital two times
	if (value_in_array(state_Q, dim, p)) {return 0;}  // creation on occupied is zero
	if (value_in_array(state_Q, dim, q)) {return 0;}  // creation on occupied is zero
	ssize_t new_state_Q[new_dim];
	new_state_Q[new_dim - 2] = p;
	new_state_Q[new_dim - 1] = q;
	for (ssize_t i = 0; i < dim; i++) {new_state_Q[i] = state_Q[i];}
	qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
	*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
	ssize_t p_pos = find_index(new_state_Q, new_dim, p);
	ssize_t q_pos = find_index(new_state_Q, new_dim, q);
	if      (p_pos < q_pos) {return n_swap_get_sign(p_pos + q_pos - 1);}
	else if (p_pos > q_pos) {return n_swap_get_sign(p_pos + q_pos);}
	else                    {return 0;} // Throwing an error
}

ssize_t smart_ca(ssize_t * P, ssize_t p, ssize_t q, ssize_t * state_Q, ssize_t dim,
                 ssize_t Q,   ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t * combo_mat)
{
	if (!value_in_array(state_Q, dim, q)) {return 0;}      // annihilate on virtual
	if (get_val_type(num_spat_orbs, num_core_elec_alpha, q) == Core)   // annihilation on frozen core, must be put back
	{
		if (p == q)  // killed a core and put it back
		{
			*P = Q;
			return 1;
		}
		else
		{return 0;} // killed a core electron, no couterpart bra to cancel this out --> ZERO
	}
	else                                                               // annihilate on occupied valence electrons
	{
		if (p == q) // killed an occupied valence and put it back
		{
			*P = Q;
			return 1;
		}
		else  // p != q
		{
			if (value_in_array(state_Q, dim, p)) {return 0;} // create on an occupied
			else                                             // create a virtual
			{
				ssize_t new_state_Q[dim];
				for(ssize_t i = 0; i < dim; i++) { if (state_Q[i] == q) {new_state_Q[i] = p;} else {new_state_Q[i] = state_Q[i];} }
				qsort(new_state_Q, dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t p_pos = find_index(new_state_Q, dim, p);
				ssize_t q_pos = find_index(state_Q,     dim, q);
				return n_swap_get_sign(p_pos + q_pos);
			}
		}
	}
}

ssize_t smart_caa(ssize_t * P,  ssize_t p, ssize_t q, ssize_t r, ssize_t * state_Q, ssize_t dim,
                  ssize_t min_num_elec, ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t * combo_mat)
{
	ssize_t n_creation = -1;
	ssize_t new_dim = dim + n_creation;
	if (new_dim < min_num_elec) {return 0;}    // too few electrons
	if (q == r) {return 0;}                    // annihilations onto the same orbital two times
	if (!value_in_array(state_Q, dim, q)) {return 0;}  // annihilate a virtual
	if (!value_in_array(state_Q, dim, r)) {return 0;}  // annihilate a virtual

	ssize_t new_state_Q[new_dim];
	ssize_t idx = 0;

	ValType p_type = get_val_type(num_spat_orbs, num_core_elec_alpha, p);  // Core electron or valence electron
	ValType q_type = get_val_type(num_spat_orbs, num_core_elec_alpha, q);
	ValType r_type = get_val_type(num_spat_orbs, num_core_elec_alpha, r);

	if (q_type == Core && r_type == Valence && p == q)  // this is  r|Q>
	{
		for (ssize_t i = 0; i < dim; i++) { if (state_Q[i] != q && state_Q[i] != r) {new_state_Q[idx++] = state_Q[i];} }
		// no need to call qsort
		*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
		ssize_t r_pos = find_index(state_Q, dim, r);
		return n_swap_get_sign(r_pos);
	}

	if (r_type == Core && q_type == Valence && p == r)  // this is -q |Q>
	{
		for (ssize_t i = 0; i < dim; i++) { if (state_Q[i] != q && state_Q[i] != r) {new_state_Q[idx++] = state_Q[i];} }
		// no need to call qsort
		*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
		ssize_t q_pos = find_index(state_Q, dim, q);
		return -1 * n_swap_get_sign(q_pos);
	}

	if (p_type == Valence && q_type == Valence && r_type == Valence)
	{
		if (value_in_array(state_Q, dim, p)) // create occupied p
		{
			if (p == q)  // now it is  r|Q>
			{
				for (ssize_t i = 0; i < dim; i++) { if (state_Q[i] != q && state_Q[i] != r){new_state_Q[idx++] = state_Q[i];} }
				new_state_Q[idx] = p;
				qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t r_pos = find_index(state_Q, dim, p);
				return n_swap_get_sign(r_pos);
			}
			else if (p == r)  // now it is  -q|Q>
			{
				for (ssize_t i = 0; i < dim; i++) { if (state_Q[i] != q && state_Q[i] != r){new_state_Q[idx++] = state_Q[i];} }
				new_state_Q[idx] = p;
				qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t q_pos = find_index(state_Q, dim, q);
				return -1 * n_swap_get_sign(q_pos);
			}
			// else --> ZERO
		}
		else  // create virtual p
		{
			for (ssize_t i = 0; i < dim; i++) { if (state_Q[i] != q && state_Q[i] != r){new_state_Q[idx++] = state_Q[i];} }
			new_state_Q[idx] = p;
			qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
			*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
			ssize_t p_pos = find_index(new_state_Q, new_dim, p);
			ssize_t q_pos = find_index(state_Q, dim, q);
			ssize_t r_pos = find_index(state_Q, dim, r);
			ssize_t n_swap = 0;
			if      (q_pos < r_pos) {n_swap = p_pos + q_pos     + r_pos;}
			else if (q_pos > r_pos) {n_swap = p_pos + q_pos - 1 + r_pos;}
			else                    {return 0;}  // ERROR
			return n_swap_get_sign(n_swap);
		}
	}
	return 0;
}

ssize_t smart_cca(ssize_t * P, ssize_t p, ssize_t q,  ssize_t r, ssize_t * state_Q, ssize_t dim,
                  ssize_t max_num_elec,   ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t * combo_mat)
{
	// 1) chekc if N+1 state exist
	ssize_t n_creation = 1;
	ssize_t new_dim = dim + n_creation;
	if (new_dim > max_num_elec) {return 0;}    // too many electrons
	if (p == q) {return 0;}                    // create two times the same orbital
	if (!value_in_array(state_Q, dim, r)) {return 0;}     // annihilate a virtual

	ssize_t new_state_Q[new_dim];
	ssize_t idx = 0;

	if (get_val_type(num_spat_orbs, num_core_elec_alpha, r) == Core)   // annihilating a core electron
	{
		if (p == r)       // This is just "-q|Q>"
	  	{
	  		if (value_in_array(state_Q, dim, q)) {return 0;} // create an occupied
			else
			{
				for (ssize_t i = 0; i < dim; i++) {new_state_Q[i] = state_Q[i];}
				new_state_Q[new_dim - 1] = q;
				qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t q_pos = find_index(new_state_Q, new_dim, q);
				return -1 * n_swap_get_sign(q_pos); // because p and q are swapped first to form "- q p r |Q>"
			}
	  	}
	  	else if (q == r)  // this is just  "p|Q>"
	  	{
	  		if (value_in_array(state_Q, dim, p)) {return 0;} // create an occupied
	  		else
	  		{
				for (ssize_t i = 0; i < dim; i++) {new_state_Q[i] = state_Q[i];}
				new_state_Q[new_dim - 1] = p;
				qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t p_pos = find_index(new_state_Q, new_dim, p);
				return n_swap_get_sign(p_pos);
			}
	  	}
	  	else
	  	{return 0;} // killed a core electron but not putting it back
	}
	else   // killing a valence electron
	{
		if (p == r) // This is just "-q|Q>"
		{
	  		if (value_in_array(state_Q, dim, q)) {return 0;} // create an occupied
	  		else
	  		{
				for (ssize_t i = 0; i < dim; i++) {new_state_Q[i] = state_Q[i];}
				new_state_Q[new_dim - 1] = q;
				qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t q_pos = find_index(new_state_Q, new_dim, q);
				return -1 * n_swap_get_sign(q_pos); // because p and q are swapped first to form "- q p r |Q>"
			}
		}
		else if (q == r) // this is just  "p|Q>"
		{
	  		if (value_in_array(state_Q, dim, p)) {return 0;} // create an occupied
	  		else
			{
				for (ssize_t i = 0; i < dim; i++) {new_state_Q[i] = state_Q[i];}
				new_state_Q[new_dim - 1] = p;
				qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t p_pos = find_index(new_state_Q, new_dim, p);
				return n_swap_get_sign(p_pos);
			}
		}
		else // this is the most complicated case...
		{
			if (value_in_array(state_Q, dim, p) || value_in_array(state_Q, dim, q)) {return 0;} // create any occupied
			else
			{
				for (ssize_t i = 0; i < dim; i++) { if (state_Q[i] != r) {new_state_Q[idx++] = state_Q[i];} }
				new_state_Q[idx++] = p;
				new_state_Q[idx]   = q;
				qsort(new_state_Q, new_dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, new_dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t p_pos = find_index(new_state_Q, new_dim, p);
				ssize_t q_pos = find_index(new_state_Q, new_dim, q);
				ssize_t r_pos = find_index(state_Q,         dim, r);
				ssize_t n_swap;
				if      (p_pos < q_pos) {n_swap = p_pos + q_pos - 1 + r_pos;}
				else if (p_pos > q_pos) {n_swap = p_pos + q_pos     + r_pos;}
				else                    {return 0;} // Throwing an error
				return n_swap_get_sign(n_swap);
			}
		}
	}
}

ssize_t smart_ccaa(ssize_t * P, ssize_t p, ssize_t q, ssize_t r, ssize_t s, ssize_t * state_Q, ssize_t dim,
                   ssize_t Q,   ssize_t num_spat_orbs, ssize_t num_core_elec_alpha, ssize_t * combo_mat)
{
	// all possiblities:
	// core      occ       virt
	// ccaa
	// ca        ca        ca
	// ca        a         c
	//           ccaa
	//           caa       c
	//           aa        cc
	if (p == q) {return 0;}  // creations     onto same orbitals
	if (r == s) {return 0;}  // annihilations onto same orbitals
	if (!value_in_array(state_Q, dim, r)) {return 0;}  // check if annihilate a virtual
	if (!value_in_array(state_Q, dim, s)) {return 0;}  // check if annihilate a virtual

	ValType p_type = get_val_type(num_spat_orbs, num_core_elec_alpha, p);   // p is core or valence electron
	ValType q_type = get_val_type(num_spat_orbs, num_core_elec_alpha, q);
	ValType r_type = get_val_type(num_spat_orbs, num_core_elec_alpha, r);
	ValType s_type = get_val_type(num_spat_orbs, num_core_elec_alpha, s);

	ssize_t new_state_Q[dim];
	ssize_t idx = 0;

	if (r_type == Core && s_type == Core)     // killing two core electrons
	{
		if (p_type == Core && q_type == Core) // all onto core electrons
		{
			if (p == r && q == s)   // one swap --> -1
			{
				*P = Q;
				return -1;
			}
			else if (p == s && q == r)  // two swaps --> 1
			{
				*P = Q;
				return 1;
			}
		}
	}

	if (r_type == Core && s_type == Valence)      // killing one core and one occupied (r is core, s is occ)
	{
		if (p == r)  // - q s |Q >
		{
			if (value_in_array(state_Q, dim, q))  // q is occ
			{
				if (q == s)
				{
					*P = Q;
					return -1;
				}
				// else -> ZERO
			}
			else // q is virt
			{
				for (ssize_t i = 0; i < dim; i++) { if(state_Q[i] != r && state_Q[i] != s) { new_state_Q[idx++] = state_Q[i]; } }
				new_state_Q[idx++] = p;
				new_state_Q[idx]   = q;
				qsort(new_state_Q, dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t q_pos = find_index(new_state_Q, dim, q);
				ssize_t s_pos = find_index(state_Q,     dim, s);
				return -1 * n_swap_get_sign(q_pos + s_pos);
			}
		}
		else if (q == r) //  p s |Q >
		{
			if (value_in_array(state_Q, dim, p))  // p is occ
			{
				if (p == s)
				{
					*P = Q;
					return 1;
				}
				// else -> ZERO
			}
			else  // p is virt
			{
				for (ssize_t i = 0; i < dim; i++) { if(state_Q[i] != r && state_Q[i] != s) { new_state_Q[idx++] = state_Q[i]; } }
				new_state_Q[idx++] = p;
				new_state_Q[idx]   = q;
				qsort(new_state_Q, dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t p_pos = find_index(new_state_Q, dim, p);
				ssize_t s_pos = find_index(state_Q,     dim, s);
				return n_swap_get_sign(p_pos + s_pos);
			}
		}
	}

	if (r_type == Valence && s_type == Core)   // killing one core and one occupied (s is core, r is occ)
	{
		if (p == s)  // q r |Q >
		{
			if (value_in_array(state_Q, dim, q))  // q is occ
			{
				if (q == r)
				{
					*P = Q;
					return 1;
				}
				// else -> ZERO
			}
			else     // q is virt
			{
				for (ssize_t i = 0; i < dim; i++) { if(state_Q[i] != r && state_Q[i] != s) { new_state_Q[idx++] = state_Q[i]; } }
				new_state_Q[idx++] = p;
				new_state_Q[idx]   = q;
				qsort(new_state_Q, dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t q_pos = find_index(new_state_Q, dim, q);
				ssize_t r_pos = find_index(state_Q,     dim, r);
				return n_swap_get_sign(q_pos + r_pos);
			}
		}
		else if (q == s)  //  - p r |Q >
		{
			if (value_in_array(state_Q, dim, p))  // p is occ
			{
				if (p == r)
				{
					*P = Q;
					return -1;
				}
				// else -> ZERO
			}
			else   // p is virt
			{
				for (ssize_t i = 0; i < dim; i++) { if(state_Q[i] != r && state_Q[i] != s) { new_state_Q[idx++] = state_Q[i]; } }
				new_state_Q[idx++] = p;
				new_state_Q[idx]   = q;
				qsort(new_state_Q, dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t p_pos = find_index(new_state_Q, dim, p);
				ssize_t r_pos = find_index(state_Q,     dim, r);
				return -1 * n_swap_get_sign(p_pos + r_pos);
			}
		}
	}

	if (r_type == Valence && s_type == Valence)  // killing two occupied orbitals
	{
		if (value_in_array(state_Q, dim, p))  // p is occ
		{
			if (value_in_array(state_Q, dim, q))  // q is occ
			{
				if (p == r && q == s)  // - p r q s |Q>
				{
					*P = Q;
					return -1;
				}
				else if (p == s && q == r) // q r p s |Q>
				{
					*P = Q;
					return 1;
				}
				// else -> ZERO
			}
			else  // q is virt
			{
				if (p == r)  // - q s |Q >
				{
					for (ssize_t i = 0; i < dim; i++) { if(state_Q[i] != r && state_Q[i] != s) { new_state_Q[idx++] = state_Q[i]; } }
					new_state_Q[idx++] = p;
					new_state_Q[idx]   = q;
					qsort(new_state_Q, dim, sizeof(ssize_t), compare);
					*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
					ssize_t q_pos = find_index(new_state_Q, dim, q);
					ssize_t s_pos = find_index(state_Q,     dim, s);
					return -1 * n_swap_get_sign(q_pos + s_pos);
				}
				else if (p == s)  // q r |Q >
				{
					for (ssize_t i = 0; i < dim; i++) { if(state_Q[i] != r && state_Q[i] != s) { new_state_Q[idx++] = state_Q[i]; } }
					new_state_Q[idx++] = p;
					new_state_Q[idx]   = q;
					qsort(new_state_Q, dim, sizeof(ssize_t), compare);
					*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
					ssize_t q_pos = find_index(new_state_Q, dim, q);
					ssize_t r_pos = find_index(state_Q,     dim, r);
					return n_swap_get_sign(q_pos + r_pos);
				}
			}

		}
		else  // p is virt
		{
			if (value_in_array(state_Q, dim, q))  // q is occ
			{
				if (q == r)  // p s |Q >
				{
					for (ssize_t i = 0; i < dim; i++) { if(state_Q[i] != r && state_Q[i] != s) { new_state_Q[idx++] = state_Q[i]; } }
					new_state_Q[idx++] = p;
					new_state_Q[idx]   = q;
					qsort(new_state_Q, dim, sizeof(ssize_t), compare);
					*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
					ssize_t p_pos = find_index(new_state_Q, dim, p);
					ssize_t s_pos = find_index(state_Q,     dim, s);
					return n_swap_get_sign(p_pos + s_pos);

				}
				else if (q == s)  // - p r |Q >
				{
					for (ssize_t i = 0; i < dim; i++) { if(state_Q[i] != r && state_Q[i] != s) { new_state_Q[idx++] = state_Q[i]; } }
					new_state_Q[idx++] = p;
					new_state_Q[idx]   = q;
					qsort(new_state_Q, dim, sizeof(ssize_t), compare);
					*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
					ssize_t p_pos = find_index(new_state_Q, dim, p);
					ssize_t r_pos = find_index(state_Q,     dim, r);
					return -1 * n_swap_get_sign(p_pos + r_pos);
				}

			}
			else  // q is virt
			{
				for (ssize_t i = 0; i < dim; i++) { if(state_Q[i] != r && state_Q[i] != s) { new_state_Q[idx++] = state_Q[i]; } }
				new_state_Q[idx++] = p;
				new_state_Q[idx]   = q;
				qsort(new_state_Q, dim, sizeof(ssize_t), compare);
				*P = find_config_index(new_state_Q, dim, num_spat_orbs, num_core_elec_alpha, combo_mat);
				ssize_t p_pos = find_index(new_state_Q, dim, p);
				ssize_t q_pos = find_index(new_state_Q, dim, q);
				ssize_t r_pos = find_index(state_Q,     dim, r);
				ssize_t s_pos = find_index(state_Q,     dim, s);
				ssize_t n_swap;
				if      (r_pos < s_pos) { n_swap = r_pos     + s_pos; }
				else if (r_pos > s_pos) { n_swap = r_pos - 1 + s_pos; }
				else                    { return 0; }
				if      (p_pos < q_pos) { n_swap += p_pos + q_pos - 1; } // q_pos comes in first and p was not there then
				else if (p_pos > q_pos) { n_swap += p_pos + q_pos;     }
				else                    { return 0; }
				return n_swap_get_sign(n_swap);
			}
		}
	}
	return 0;
}

/* end of "b" operator computations */



/* Interfaces for Python ctypes.cdll.LoadLibraray(...) */

void build_A(ssize_t * configlist[], double * zmatlist[],          double  CCAA_IJ_storage[],
             ssize_t I,              ssize_t J,                    ssize_t num_charge_state,
             ssize_t num_spat_orbs,  ssize_t num_core_elec_alpha,  ssize_t * combo_mat[],
             ssize_t num_elec[],     ssize_t num_states[],         ssize_t num_states_used[],  ssize_t nthd)
{
	ssize_t P;
	ssize_t expect_val;
	ssize_t num_spin_orb = 2 * num_spat_orbs;
	ssize_t min_num_elec = find_min(num_elec, num_charge_state);
	ssize_t tensor_size  = num_spin_orb;
	ssize_t dim = num_elec[J];

	if (num_elec[I] - num_elec[J] == -1)  // Just a guard
	{
		for (ssize_t Q = 0; Q < num_states[J]; Q++)
		{
			for (ssize_t p = 0; p < num_spin_orb; p++)
			{
				expect_val = smart_a(&P, p, &( configlist[J][dim * Q] ), dim, min_num_elec, num_spat_orbs, num_core_elec_alpha, combo_mat[I]);
				if (expect_val != 0) // non-zero case
				{
					for (ssize_t i = 0; i < num_states_used[I]; i++)
					{
						for (ssize_t j = 0; j < num_states_used[J]; j++)
						{
							CCAA_IJ_storage[ (i*num_states_used[J] + j)*tensor_size + p ] +=
							                zmatlist[I][num_states[I] * i + P] * zmatlist[J][num_states[J] * j + Q] * expect_val;
						}
					}
				}
			}
		}
	}
}

void build_C(ssize_t * configlist[], double * zmatlist[],          double  CCAA_IJ_storage[],
             ssize_t I,              ssize_t J,                    ssize_t num_charge_state,
             ssize_t num_spat_orbs,  ssize_t num_core_elec_alpha,  ssize_t * combo_mat[],
             ssize_t num_elec[],     ssize_t num_states[],         ssize_t num_states_used[],  ssize_t nthd)
{
	ssize_t P;
	ssize_t expect_val;
	ssize_t num_spin_orb = 2 * num_spat_orbs;
	ssize_t max_num_elec = find_max(num_elec, num_charge_state);
	ssize_t tensor_size  = num_spin_orb;
	ssize_t dim = num_elec[J];

	if (num_elec[I] - num_elec[J] == 1)  // Just a guard
	{
		for (ssize_t Q = 0; Q < num_states[J]; Q++)
		{
			for (ssize_t p = 0; p < num_spin_orb; p++)
			{
				expect_val = smart_c(&P, p, &( configlist[J][dim * Q] ), dim, max_num_elec, num_spat_orbs, num_core_elec_alpha, combo_mat[I]);
				if (expect_val != 0) // non-zero case
				{
					for (ssize_t i = 0; i < num_states_used[I]; i++)
					{
						for (ssize_t j = 0; j < num_states_used[J]; j++)
						{
							CCAA_IJ_storage[ (i*num_states_used[J] + j)*tensor_size + p ] +=
										zmatlist[I][num_states[I] * i + P] * zmatlist[J][num_states[J] * j + Q] * expect_val;
						}
					}
				}
			}
		}
	}
}

void build_AA(ssize_t * configlist[], double * zmatlist[],          double  CCAA_IJ_storage[],
              ssize_t I,              ssize_t J,                    ssize_t num_charge_state,
              ssize_t num_spat_orbs,  ssize_t num_core_elec_alpha,  ssize_t * combo_mat[],
              ssize_t num_elec[],     ssize_t num_states[],         ssize_t num_states_used[],  ssize_t nthd)
{
	ssize_t P;
	ssize_t expect_val;
	ssize_t num_spin_orb = 2 * num_spat_orbs;
	ssize_t tensor_size  = num_spin_orb * num_spin_orb;
	ssize_t min_num_elec = find_min(num_elec, num_charge_state);
	ssize_t dim = num_elec[J];

	if (num_elec[I] - num_elec[J] == -2)  // Just a guard
	{
		for (ssize_t Q = 0; Q < num_states[J]; Q++)
		{
			for (ssize_t p = 0; p < num_spin_orb; p++)
			{
				for (ssize_t q = 0; q < num_spin_orb; q++)
				{
					expect_val = smart_aa(&P, p, q, &( configlist[J][dim * Q] ), dim, min_num_elec, num_spat_orbs, num_core_elec_alpha, combo_mat[I]);
					if (expect_val != 0) // non-zero case
					{
						for (ssize_t i = 0; i < num_states_used[I]; i++)
						{
							for (ssize_t j = 0; j < num_states_used[J]; j++)
							{
								CCAA_IJ_storage[ (i*num_states_used[J] + j)*tensor_size + num_spin_orb * p + q ] +=
										zmatlist[I][num_states[I] * i + P] * zmatlist[J][num_states[J] * j + Q] * expect_val;
							}
						}
					}
				}
			}
		}
	}
}

void build_CC(ssize_t * configlist[], double * zmatlist[],          double  CCAA_IJ_storage[],
              ssize_t I,              ssize_t J,                    ssize_t num_charge_state,
              ssize_t num_spat_orbs,  ssize_t num_core_elec_alpha,  ssize_t * combo_mat[],
              ssize_t num_elec[],     ssize_t num_states[],         ssize_t num_states_used[],  ssize_t nthd)
{
	ssize_t P;
	ssize_t expect_val;
	ssize_t num_spin_orb = 2 * num_spat_orbs;
	ssize_t tensor_size  = num_spin_orb * num_spin_orb;
	ssize_t max_num_elec = find_max(num_elec, num_charge_state);
	ssize_t dim = num_elec[J];

	if (num_elec[I] - num_elec[J] == 2)  // Just a guard
	{
		for (ssize_t Q = 0; Q < num_states[J]; Q++)
		{
			for (ssize_t p = 0; p < num_spin_orb; p++)
			{
				for (ssize_t q = 0; q < num_spin_orb; q++)
				{
					expect_val = smart_cc(&P, p, q, &( configlist[J][dim * Q] ), dim, max_num_elec, num_spat_orbs, num_core_elec_alpha, combo_mat[I]);
					if (expect_val != 0) // non-zero case
					{
						for (ssize_t i = 0; i < num_states_used[I]; i++)
						{
							for (ssize_t j = 0; j < num_states_used[J]; j++)
							{
								CCAA_IJ_storage[ (i*num_states_used[J] + j)*tensor_size + num_spin_orb * p + q ] += 
											zmatlist[I][num_states[I] * i + P] * zmatlist[J][num_states[J] * j + Q] * expect_val;
							}
						}
					}
				}
			}
		}
	}

}

void build_CA(ssize_t * configlist[], double * zmatlist[],         double  CCAA_IJ_storage[],
               ssize_t I,              ssize_t J,
               ssize_t num_spat_orbs,  ssize_t num_core_elec_alpha, ssize_t * combo_mat[],
               ssize_t num_elec[],     ssize_t num_states[],        ssize_t num_states_used[],  ssize_t nthd)
{
	ssize_t P;
	ssize_t expect_val;
	ssize_t num_spin_orb = 2 * num_spat_orbs;
	ssize_t tensor_size  = num_spin_orb * num_spin_orb;
	ssize_t dim = num_elec[J];

	if (num_elec[I] - num_elec[J] == 0)  // Just a guard
	{
		for (ssize_t Q = 0; Q < num_states[J]; Q++)
		{
			for (ssize_t p = 0; p < num_spin_orb; p++)
			{
				for (ssize_t q = 0; q < num_spin_orb; q++)
				{
					expect_val = smart_ca(&P, p, q, &( configlist[J][dim * Q] ), dim, Q, num_spat_orbs, num_core_elec_alpha, combo_mat[I]);
					if (expect_val != 0) // non-zero case
					{
						for (ssize_t i = 0; i < num_states_used[I]; i++)
						{
							for (ssize_t j = 0; j < num_states_used[J]; j++)
							{
								CCAA_IJ_storage[ (i*num_states_used[J] + j)*tensor_size + num_spin_orb * p + q ] += 
											zmatlist[I][num_states[I] * i + P] * zmatlist[J][num_states[J] * j + Q] * expect_val;
							}
						}
					}
				}
			}
		}
	}
}

void build_CAA(ssize_t * configlist[], double * zmatlist[],          double  CCAA_IJ_storage[],
               ssize_t I,              ssize_t J,                    ssize_t num_charge_state,
               ssize_t num_spat_orbs,  ssize_t num_core_elec_alpha,  ssize_t * combo_mat[],
               ssize_t num_elec[],     ssize_t num_states[],         ssize_t num_states_used[],  ssize_t nthd)
{
	ssize_t P;
	ssize_t expect_val;
	ssize_t num_spin_orb = 2 * num_spat_orbs;
	ssize_t tensor_size  = num_spin_orb * num_spin_orb * num_spin_orb;
	ssize_t min_num_elec = find_min(num_elec, num_charge_state);
	ssize_t dim = num_elec[J];

	if (num_elec[I] - num_elec[J] == -1)  // Just a guard
	{
		for (ssize_t Q = 0; Q < num_states[J]; Q++)
		{
			for (ssize_t p = 0; p < num_spin_orb; p++)
			{
				for (ssize_t q = 0; q < num_spin_orb; q++)
				{
					for (ssize_t r = 0; r < num_spin_orb; r++)
					{
						expect_val = smart_caa(&P, p, q, r, &( configlist[J][dim * Q] ), dim, min_num_elec, num_spat_orbs, num_core_elec_alpha, combo_mat[I]);
						if (expect_val != 0) // non-zero case
						{
							for (ssize_t i = 0; i < num_states_used[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used[J]; j++)
								{
									CCAA_IJ_storage[ (i*num_states_used[J] + j)*tensor_size + num_spin_orb * (num_spin_orb * p + q) + r ] += 
												zmatlist[I][num_states[I] * i + P] * zmatlist[J][num_states[J] * j + Q] * expect_val;
								}
							}
						}
					}
				}
			}
		}
	}
}

void build_CCA(ssize_t * configlist[], double * zmatlist[],          double  CCAA_IJ_storage[],
               ssize_t I,              ssize_t J,                    ssize_t num_charge_state,
               ssize_t num_spat_orbs,  ssize_t num_core_elec_alpha,  ssize_t * combo_mat[],
               ssize_t num_elec[],     ssize_t num_states[],         ssize_t num_states_used[],  ssize_t nthd)
{
	ssize_t P;
	ssize_t expect_val;
	ssize_t num_spin_orb = 2 * num_spat_orbs;
	ssize_t max_num_elec = find_max(num_elec, num_charge_state);
	ssize_t tensor_size  = num_spin_orb * num_spin_orb * num_spin_orb;
	ssize_t dim = num_elec[J];

	if (num_elec[I] - num_elec[J] == 1)  // Just a guard
	{
		for (ssize_t Q = 0; Q < num_states[J]; Q++)
		{
			for (ssize_t p = 0; p < num_spin_orb; p++)
			{
				for (ssize_t q = 0; q < num_spin_orb; q++)
				{
					for (ssize_t r = 0; r < num_spin_orb; r++)
					{
						expect_val = smart_cca(&P, p, q, r, &( configlist[J][dim * Q] ), dim, max_num_elec, num_spat_orbs, num_core_elec_alpha, combo_mat[I]);
						if (expect_val != 0) // non-zero case
						{
							for (ssize_t i = 0; i < num_states_used[I]; i++)
							{
								for (ssize_t j = 0; j < num_states_used[J]; j++)
								{
									CCAA_IJ_storage[ (i*num_states_used[J] + j)*tensor_size + num_spin_orb * (num_spin_orb * p + q) + r ] += 
												zmatlist[I][num_states[I] * i + P] * zmatlist[J][num_states[J] * j + Q] * expect_val;
								}
							}
						}
					}
				}
			}
		}
	}
}

void build_CCAA(ssize_t * configlist[], double * zmatlist[],          double  CCAA_IJ_storage[],
                ssize_t I,              ssize_t J,
                ssize_t num_spat_orbs,  ssize_t num_core_elec_alpha,  ssize_t * combo_mat[],
                ssize_t num_elec[],     ssize_t num_states[],         ssize_t num_states_used[],  ssize_t nthd)
{
	ssize_t P;
	ssize_t expect_val;
	ssize_t num_spin_orb = 2 * num_spat_orbs;
	ssize_t tensor_size  = num_spin_orb * num_spin_orb * num_spin_orb * num_spin_orb;
	ssize_t dim = num_elec[J];

	if (num_elec[I] - num_elec[J] == 0)  // Just a guard
	{
		for (ssize_t Q = 0; Q < num_states[J]; Q++)
		{
			for (ssize_t p = 0; p < num_spin_orb; p++)
			{
				for (ssize_t q = 0; q < num_spin_orb; q++)
				{
					for (ssize_t r = 0; r < num_spin_orb; r++)
					{
						for (ssize_t s = 0; s < num_spin_orb; s++)
						{
							expect_val = smart_ccaa(&P, p, q, r, s, &( configlist[J][dim * Q] ), dim, Q, num_spat_orbs, num_core_elec_alpha, combo_mat[I]);
							if (expect_val != 0) // non-zero case
							{
								for (ssize_t i = 0; i < num_states_used[I]; i++)
								{
									for (ssize_t j = 0; j < num_states_used[J]; j++)
									{
										CCAA_IJ_storage[ (i*num_states_used[J] + j)*tensor_size + num_spin_orb * (num_spin_orb * (num_spin_orb * p + q) + r) + s ] +=
														zmatlist[I][num_states[I] * i + P] * zmatlist[J][num_states[J] * j + Q] * expect_val;
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
