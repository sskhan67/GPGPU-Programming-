import numpy

def one_electron(O_raw):
    n_orb = O_raw.shape[0]
    O = numpy.zeros((2*n_orb, 2*n_orb))
    for p in range(n_orb):
        for q in range(n_orb):
            O[p, q]             = O_raw[p, q]
            O[p+n_orb, q+n_orb] = O_raw[p, q]
    return O

def two_electron(O_raw):
    n_orb = O_raw.shape[0]
    O = numpy.zeros((2*n_orb, 2*n_orb, 2*n_orb, 2*n_orb))
    for p in range(n_orb):
        for q in range(n_orb):
            for r in range(n_orb):
                for s in range(n_orb):
                    O[p, q, r, s]                         = O_raw[p, q, r, s]
                    O[p+n_orb, q, r+n_orb, s]             = O_raw[p, q, r, s]
                    O[p, q+n_orb, r, s+n_orb]             = O_raw[p, q, r, s]
                    O[p+n_orb, q+n_orb, r+n_orb, s+n_orb] = O_raw[p, q, r, s]
    return O
