import numpy
import psi4

def brabraketket(braketbraket):
    n_orb = braketbraket.shape[0]
    V = numpy.zeros((n_orb, n_orb, n_orb, n_orb))
    for p in range(n_orb):
        for q in range(n_orb):
            for r in range(n_orb):
                for s in range(n_orb):
                    V[p,q,r,s] = braketbraket[p,r,q,s]
    return V

class supplementary(object):
    def __init__(self, mol, bas, ints):
        self.mol, self.bas, self.ints = mol, bas, ints

def AO_ints(z_matrix, basis, max_mem=1e9):
    mol  = psi4.geometry(z_matrix)
    mol.update_geometry()
    bas  = psi4.core.BasisSet.build(mol, target=basis)
    ints = psi4.core.MintsHelper(bas)
    mem = ints.nbf()**4 * 8	# assuming 8 bytes per integral
    if mem>max_mem:  raise ValueError("Raise max_mem to above {} (bytes) to store this many integrals.".format(mem))
    S = numpy.array(ints.ao_overlap())
    U = numpy.array(ints.ao_potential())
    T = numpy.array(ints.ao_kinetic())
    V = brabraketket(numpy.array(ints.ao_eri()))
    return S,T,U,V,supplementary(mol,bas,ints)
