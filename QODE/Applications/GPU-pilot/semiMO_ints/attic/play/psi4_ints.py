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

class AO_ints(object):
    def __init__(self, z_matrix, basis, max_mem=1e9):
        mol  = psi4.geometry(z_matrix)
        mol.update_geometry()
        bas  = psi4.core.BasisSet.build(mol, target=basis)
        self.ints = psi4.core.MintsHelper(bas)
    def overlap(self):
        return numpy.array(self.ints.ao_overlap())
    def potential(self):
        return numpy.array(self.ints.ao_potential())
    def kinetic(self):
        return numpy.array(self.ints.ao_kinetic())
    def electron_repulsion(self):
        mem = self.ints.nbf()**4 * 8	# assuming 8 bytes per integral
        if mem>max_mem:  raise ValueError("Raise max_mem to above {} (bytes) to store this many integrals.".format(mem))
        return brabraketket(numpy.array(self.ints.ao_eri()))
    def info(self):
        return supplementary(mol,bas,ints)
    @staticmethod
    def all(z_matrix, basis, max_mem=1e9):
        integrals = AO_ints(z_matrix, basis, max_mem=1e9)
        return integrals.overlap(), integrals.potential(), integrals.kinetic(), integrals.electron_repulsion(), integrals.info()
