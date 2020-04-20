import numpy as np
# this file contains all the intregral transformations performed 


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

def get_semiMO_basis(S, C):
    Be_coef = len(C)
    n_Be    = int(len(S)/Be_coef)
    l = 0
    for n in range(n_Be):
        beg   = int(n*Be_coef)
        end   = int(beg + Be_coef)
        S_Be  = S[:,beg:end]
        Spp_slice   = S_Be @ C
        if l == 0:
            Spp = Spp_slice
            l  += 1
        else:
            Spp = np.concatenate( (Spp, Spp_slice), axis = 1)
    m = 0
    for n in range(n_Be):
        beg   = int(n*Be_coef)
        end   = int(beg+ Be_coef)
        S_Be  = Spp[beg:end,:]
        S_semiMO_slice = C.T @ S_Be
        if m == 0:
            S_semiMO = S_semiMO_slice
            m += 1
        else:
            S_semiMO = np.concatenate( (S_semiMO, S_semiMO_slice), axis = 0)
    return(S_semiMO)


def get_G_semiMO_basis(g_c, Be_C):
    for dimension in range(4):
        l = 0
        Be_coef = len(Be_C)
        n_Be    = int(len(g_c)/Be_coef)
        for n in range(n_Be):
            beg = int(n*Be_coef)
            end = int(beg + Be_coef)
            if   dimension == 0:
                #print('beg=', beg, 'end=', end)
                g_c_Be   = g_c[beg:end,:,:,:]
            elif dimension == 1:
                g_c_Be   = G0[:,beg:end,:,:]
            elif dimension == 2:
                g_c_Be   = G1[:,:,beg:end,:]
            elif dimension == 3:
                g_c_Be   = G2[:,:,:,beg:end]
            g_c_slice   = np.tensordot(Be_C, g_c_Be, axes=([0],[dimension]))
            g_c_slice   = np.moveaxis(g_c_slice, 0, dimension)
            if l == 0:
                if dimension == 0:
                    G0 = g_c_slice
                if dimension == 1:
                    G1 = g_c_slice
                if dimension == 2:
                    G2 = g_c_slice
                if dimension == 3:
                    G3 = g_c_slice
                #l += 1
            else:
                if dimension == 0:
                    G0 = np.concatenate( (G0, g_c_slice), axis = dimension)
                if dimension == 1:
                    G1 = np.concatenate( (G1, g_c_slice), axis = dimension)
                if dimension == 2:
                    G2 = np.concatenate( (G2, g_c_slice), axis = dimension)
                if dimension == 3:
                    G3 = np.concatenate( (G3, g_c_slice), axis = dimension)
            l += 1
    return(G3)



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

