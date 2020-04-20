import numpy as np

def alpha(index):
    return index*2

def beta(index):
    return index*2 + 1 

def H_spatial_to_spin_orb_grid(H):
    n_spatial  = len(H)
    n_spin_orb = len(H)*2
    Hcore      = np.zeros((n_spin_orb, n_spin_orb))
    for p in range(n_spatial):
        for q in range(n_spatial):
            Hcore[alpha(p), alpha(q)] = H[p, q]
            Hcore[beta(p) , beta(q) ] = H[p, q]
    return Hcore

def g_spatial_to_spin_orb_chem_grid(g):
    n_spatial  = len(g)
    n_spin_orb = n_spatial*2
    v          = np.zeros((n_spin_orb, n_spin_orb, n_spin_orb, n_spin_orb))
    for p in range(n_spatial):
        for q in range(n_spatial):
            for r in range(n_spatial):
                for s in range(n_spatial):
                    v[alpha(p),alpha(q),alpha(r),alpha(s)]  = g[p,q,r,s]
                    v[alpha(p),alpha(q), beta(r), beta(s)]  = g[p,q,r,s]
                    v[beta( p), beta(q),alpha(r),alpha(s)]  = g[p,q,r,s]
                    v[beta( p), beta(q), beta(r), beta(s)]  = g[p,q,r,s]
    return v

def g_spatial_to_spin_orb_phys_grid(g):
    n_spatial  = len(g)
    n_spin_orb = n_spatial*2
    v          = np.zeros((n_spin_orb, n_spin_orb, n_spin_orb, n_spin_orb))
    for p in range(n_spatial):
        for q in range(n_spatial):
            for r in range(n_spatial):
                for s in range(n_spatial):
                    v[alpha(p),alpha(r),alpha(q),alpha(s)]  = g[p,q,r,s]
                    v[alpha(p), beta(r),alpha(q), beta(s)]  = g[p,q,r,s]
                    v[beta( p),alpha(r), beta(q),alpha(s)]  = g[p,q,r,s]
                    v[beta( p), beta(r), beta(q), beta(s)]  = g[p,q,r,s]
    return v



def one_electron_to_spin_block(O_raw):
    n_orb = O_raw.shape[0]
    O = np.zeros((2*n_orb, 2*n_orb))
    for p in range(n_orb):
        for q in range(n_orb):
            O[p, q]             = O_raw[p, q]
            O[p+n_orb, q+n_orb] = O_raw[p, q]
    return O


def two_electron_to_spin_phys_block(O_raw):
    n_orb = O_raw.shape[0]
    O = np.zeros((2*n_orb, 2*n_orb, 2*n_orb, 2*n_orb))
    for p in range(n_orb):
        for q in range(n_orb):
            for r in range(n_orb):
                for s in range(n_orb):
                    O[p, q, r, s]                         = O_raw[p, r, q, s]
                    O[p+n_orb, q, r+n_orb, s]             = O_raw[p, r, q, s]
                    O[p, q+n_orb, r, s+n_orb]             = O_raw[p, r, q, s]
                    O[p+n_orb, q+n_orb, r+n_orb, s+n_orb] = O_raw[p, r, q, s]
    return O


def two_electron_to_spin_chem_block(O_raw):
    n_orb = O_raw.shape[0]
    O = np.zeros((2*n_orb, 2*n_orb, 2*n_orb, 2*n_orb))
    for p in range(n_orb):
        for q in range(n_orb):
            for r in range(n_orb):
                for s in range(n_orb):
                    O[p,       r,       q,       s]       = O_raw[p, r, q, s]
                    O[p+n_orb, r+n_orb, q,       s]       = O_raw[p, r, q, s]
                    O[p,       r,       q+n_orb, s+n_orb] = O_raw[p, r, q, s]
                    O[p+n_orb, r+n_orb, q+n_orb, s+n_orb] = O_raw[p, r, q, s]
    return O



'''
def g_spatial_to_spin_orb_phys(g):
    n_spatial  = len(g)
    n_spin_orb = n_spatial*2
    v          = np.zeros((n_spin_orb, n_spin_orb, n_spin_orb, n_spin_orb))
    for p in range(n_spatial):
        for q in range(n_spatial):
            for r in range(n_spatial):
                for s in range(n_spatial):
                    v[alpha(p),alpha(r),alpha(q),alpha(s)]  = g[p,q,r,s]
                    v[alpha(p),alpha(r),beta( q),beta( s)]  = g[p,q,r,s]
                    v[beta( p),beta( r),alpha(q),alpha(s)]  = g[p,q,r,s]
                    v[beta( p),beta( r),beta( q),beta( s)]  = g[p,q,r,s]
    return v
'''

