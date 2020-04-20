import numpy
import Double


class empty(object):  pass

def monomers(sys_list):
    ints1 = []
    for system in sys_list:
        S, T, U, V = integrals(system)
        ints = empty()
        ints.S = S
        ints.T = T
        ints.U = U
        ints.V = V
        ints1 += [ints]
    return ints1

def dimers(sys_list, ints1):
    ints2 = []
    for nA,sysA in enumerate(sys_list):
        intsA = ints1[nA]
        dimA = intsA.S.shape[0]
        intsAb = []
        for nB,sysB in enumerate(sys_list):
            intsB = ints1[nB]
            if nB<=nA:
                intsAb += [None]
            else:
                system = combine(sysA,sysB)
                Sraw, Traw, Uraw, Vraw = integrals(system)
                intsAB = empty()
                intsAB.S = [[None,None],[None,None]]
                intsAB.S[0][1] = numpy.array(Sraw[:dimA,dimA:])
                intsAB.S[1][0] = numpy.array(Sraw[dimA:,:dimA])
                intsAB.T = [[None,None],[None,None]]
                intsAB.T[0][1] = numpy.array(Traw[:dimA,dimA:])
                intsAB.T[1][0] = numpy.array(Traw[dimA:,:dimA])
                intsAB.U = [[[None,None],[None,None]],[[None,None],[None,None]]]
                intsAB.U[0][0][1] = Uraw[:dimA,dimA:] / 2
                intsAB.U[0][1][0] = Uraw[dimA:,:dimA] / 2
                intsAB.U[0][1][1] = Uraw[dim:,dim:] - intsB.U
                intsAB.U[1][0][0] = Uraw[:dim,:dim] - intsA.U
                intsAB.U[1][0][1] = Uraw[:dimA,dimA:] / 2
                intsAB.U[1][1][0] = Uraw[dimA:,:dimA] / 2
                intsAB.V = [[[[None,None],[None,None]],[[None,None],[None,None]]],[[[None,None],[None,None]],[[None,None],[None,None]]]]
                intsAB.V[0][0][0][1] = numpy.array(Vraw[:dimA,:dimA,:dimA,dimA:])
                intsAB.V[0][0][1][0] = numpy.array(Vraw[:dimA,:dimA,dimA:,:dimA])
                intsAB.V[0][0][1][1] = numpy.array(Vraw[:dimA,:dimA,dimA:,dimA:])
                intsAB.V[0][1][0][0] = numpy.array(Vraw[:dimA,dimA:,:dimA,:dimA])
                intsAB.V[0][1][0][1] = numpy.array(Vraw[:dimA,dimA:,:dimA,dimA:])
                intsAB.V[0][1][1][0] = numpy.array(Vraw[:dimA,dimA:,dimA:,:dimA])
                intsAB.V[0][1][1][1] = numpy.array(Vraw[:dimA,dimA:,dimA:,dimA:])
                intsAB.V[1][0][0][0] = numpy.array(Vraw[dimA:,:dimA,:dimA,:dimA])
                intsAB.V[1][0][0][1] = numpy.array(Vraw[dimA:,:dimA,:dimA,dimA:])
                intsAB.V[1][0][1][0] = numpy.array(Vraw[dimA:,:dimA,dimA:,:dimA])
                intsAB.V[1][0][1][1] = numpy.array(Vraw[dimA:,:dimA,dimA:,dimA:])
                intsAB.V[1][1][0][0] = numpy.array(Vraw[dimA:,dimA:,:dimA,:dimA])
                intsAB.V[1][1][0][1] = numpy.array(Vraw[dimA:,dimA:,:dimA,dimA:])
                intsAB.V[1][1][1][0] = numpy.array(Vraw[dimA:,dimA:,dimA:,:dimA])
                intsAb += [intsAB]
        ints2 += [intsAb]
    return ints2

def trimers(sys_list, ints1, ints2):
    ints3 = []
    for nA,sysA in enumerate(sys_list):
        dimA = ints1[nA].S.shape[0]
        intsAbc = []
        for nB,sysB in enumerate(sys_list):
            dimAB = dimA + ints1[nB].S.shape[0]
            intsAB = ints2[nA][nB]
            intsABc = []
            for nC,sysC in enumerate(sys_list):
                intsAC = ints2[nA][nC]
                intsBC = ints2[nB][nC]
                if nC<=nB or nB<=nA:
                    intsABc += [None]
                else:
                    system = combine(sysA,sysB,sysC)
                    Sraw, Traw, Uraw, Vraw = integrals(system)
                    intsABC = empty()
                    intsABC.U = [[[None,None,None],[None,None,None],[None,None,None]],[[None,None,None],[None,None,None],[None,None,None]],[[None,None,None],[None,None,None],[None,None,None]]]
                    intsABC.U[0][1][2] = Uraw[dimA:dimAB,dimAB:] - intsBC.U[0][0][1] - intsBC.U[1][0][1]
                    intsABC.U[0][2][1] = Uraw[dimAB:,dimA:dimAB] - intsBC.U[0][1][0] - intsBC.U[1][1][0]
                    intsABC.U[1][0][2] = Uraw[:dimA,dimAB:]      - intsAC.U[0][0][1] - intsAC.U[1][0][1]
                    intsABC.U[1][2][0] = Uraw[dimAB:,:dimA]      - intsAC.U[0][1][0] - intsAC.U[1][1][0]
                    intsABC.U[2][0][1] = Uraw[:dimA,dimA:dimAB]  - intsAB.U[0][0][1] - intsAB.U[1][0][1]
                    intsABC.U[2][1][0] = Uraw[dimA:dimAB,:dimA]  - intsAB.U[0][1][0] - intsAB.U[1][1][0]
                    intsABC.V = [[[[None,None,None],[None,None,None],[None,None,None]],[[None,None,None],[None,None,None],[None,None,None]],[[None,None,None],[None,None,None],[None,None,None]]],[[[None,None,None],[None,None,None],[None,None,None]],[[None,None,None],[None,None,None],[None,None,None]],[[None,None,None],[None,None,None],[None,None,None]]],[[[None,None,None],[None,None,None],[None,None,None]],[[None,None,None],[None,None,None],[None,None,None]],[[None,None,None],[None,None,None],[None,None,None]]]]
                    intsABC.V[0][0][1][2] = numpy.array(Vraw[:dimA,:dimA,dimA:dimAB,dimAB:])
                    intsABC.V[0][0][2][1] = numpy.array(Vraw[:dimA,:dimA,dimAB:,dimA:dimAB])
                    intsABC.V[0][1][0][2] = numpy.array(Vraw[:dimA,dimA:dimAB,:dimA,dimAB:])
                    intsABC.V[0][2][0][1] = numpy.array(Vraw[:dimA,dimAB:,:dimA,dimA:dimAB])
                    intsABC.V[0][1][2][0] = numpy.array(Vraw[:dimA,dimA:dimAB,dimAB:,:dimA])
                    intsABC.V[0][2][1][0] = numpy.array(Vraw[:dimA,dimAB:,dimA:dimAB,:dimA])
                    intsABC.V[1][0][0][2] = numpy.array(Vraw[dimA:dimAB,:dimA,:dimA,dimAB:])
                    intsABC.V[2][0][0][1] = numpy.array(Vraw[dimAB:,:dimA,:dimA,dimA:dimAB])
                    intsABC.V[1][0][2][0] = numpy.array(Vraw[dimA:dimAB,:dimA,dimAB:,:dimA])
                    intsABC.V[2][0][1][0] = numpy.array(Vraw[dimAB:,:dimA,dimA:dimAB,:dimA])
                    intsABC.V[1][2][0][0] = numpy.array(Vraw[dimA:dimAB,dimAB:,:dimA,:dimA])
                    intsABC.V[2][1][0][0] = numpy.array(Vraw[dimAB:,dimA:dimAB,:dimA,:dimA])
                    intsABC.V[1][1][0][2] = numpy.array(Vraw[dimA:dimAB,dimA:dimAB,:dimA,dimAB:])
                    intsABC.V[1][1][2][0] = numpy.array(Vraw[dimA:dimAB,dimA:dimAB,dimAB:,:dimA])
                    intsABC.V[1][0][1][2] = numpy.array(Vraw[dimA:dimAB,:dimA,dimA:dimAB,dimAB:])
                    intsABC.V[1][2][1][0] = numpy.array(Vraw[dimA:dimAB,dimAB:,dimA:dimAB,:dimA])
                    intsABC.V[1][0][2][1] = numpy.array(Vraw[dimA:dimAB,:dimA,dimAB:,dimA:dimAB])
                    intsABC.V[1][2][0][1] = numpy.array(Vraw[dimA:dimAB,dimAB:,:dimA,dimA:dimAB])
                    intsABC.V[0][1][1][2] = numpy.array(Vraw[:dimA,dimA:dimAB,dimA:dimAB,dimAB:])
                    intsABC.V[2][1][1][0] = numpy.array(Vraw[dimAB:,dimA:dimAB,dimA:dimAB,:dimA])
                    intsABC.V[0][1][2][1] = numpy.array(Vraw[:dimA,dimA:dimAB,dimAB:,dimA:dimAB])
                    intsABC.V[2][1][0][1] = numpy.array(Vraw[dimAB:,dimA:dimAB,:dimA,dimA:dimAB])
                    intsABC.V[0][2][1][1] = numpy.array(Vraw[:dimA,dimAB:,dimA:dimAB,dimA:dimAB])
                    intsABC.V[2][0][1][1] = numpy.array(Vraw[dimAB:,:dimA,dimA:dimAB,dimA:dimAB])
                    intsABC.V[2][2][0][1] = numpy.array(Vraw[dimAB:,dimAB:,:dimA,dimA:dimAB])
                    intsABC.V[2][2][1][0] = numpy.array(Vraw[dimAB:,dimAB:,dimA:dimAB,:dimA])
                    intsABC.V[2][0][2][1] = numpy.array(Vraw[dimAB:,:dimA,dimAB:,dimA:dimAB])
                    intsABC.V[2][1][2][0] = numpy.array(Vraw[dimAB:,dimA:dimAB,dimAB:,:dimA])
                    intsABC.V[2][0][1][2] = numpy.array(Vraw[dimAB:,:dimA,dimA:dimAB,dimAB:])
                    intsABC.V[2][1][0][2] = numpy.array(Vraw[dimAB:,dimA:dimAB,:dimA,dimAB:])
                    intsABC.V[0][2][2][1] = numpy.array(Vraw[:dimA,dimAB:,dimAB:,dimA:dimAB])
                    intsABC.V[1][2][2][0] = numpy.array(Vraw[dimA:dimAB,dimAB:,dimAB:,:dimA])
                    intsABC.V[0][2][1][2] = numpy.array(Vraw[:dimA,dimAB:,dimA:dimAB,dimAB:])
                    intsABC.V[1][2][0][2] = numpy.array(Vraw[dimA:dimAB,dimAB:,:dimA,dimAB:])
                    intsABC.V[0][1][2][2] = numpy.array(Vraw[:dimA,dimA:dimAB,dimAB:,dimAB:])
                    intsABC.V[1][0][2][2] = numpy.array(Vraw[dimA:dimAB,:dimA,dimAB:,dimAB:])
                    intsABc += [intsABC]
            intsAbc += [intsABc]
        ints3 += [intsAbc]
    return ints3



def tetramers(sys_list, ints1):
    ints4 = []
    for nA,sysA in enumerate(sys_list):
        dimA = ints1[nA].S.shape[0]
        intsAbcd = []
        for nB,sysB in enumerate(sys_list):
            dimAB = dimA + ints1[nB].S.shape[0]
            intsABcd = []
            for nC,sysC in enumerate(sys_list):
                dimABC = dimAB + ints1[nC].S.shape[0]
                intsABCd = []
                for nD,sysD in enumerate(sys_list):
                    if nD<=nC or nC<=nB or nB<=nA:
                        intsABCd += [None]
                    else:
                        system = combine(sysA,sysB,sysC,sysD)
                        Sraw, Traw, Uraw, Vraw = integrals(system)
                        intsABCD = empty()
                        intsABCD.V = [[[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]]],[[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]]],[[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]]],[[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]],[[None,None,None,None],[None,None,None,None],[None,None,None,None],[None,None,None,None]]]]
                        intsABCD.V[0][1][2][3] = numpy.array(Vraw[:dimA,dimA:dimAB,dimAB:dimABC,dimABC:])
                        intsABCD.V[0][1][3][2] = numpy.array(Vraw[:dimA,dimA:dimAB,dimABC:,dimAB:dimABC])
                        intsABCD.V[0][2][1][3] = numpy.array(Vraw[:dimA,dimAB:dimABC,dimA:dimAB,dimABC:])
                        intsABCD.V[0][3][1][2] = numpy.array(Vraw[:dimA,dimABC:,dimA:dimAB,dimAB:dimABC])
                        intsABCD.V[0][2][3][1] = numpy.array(Vraw[:dimA,dimAB:dimABC,dimABC:,dimA:dimAB])
                        intsABCD.V[0][3][2][1] = numpy.array(Vraw[:dimA,dimABC:,dimAB:dimABC,dimA:dimAB])
                        intsABCD.V[1][0][2][3] = numpy.array(Vraw[dimA:dimAB,:dimA,dimAB:dimABC,dimABC:])
                        intsABCD.V[1][0][3][2] = numpy.array(Vraw[dimA:dimAB,:dimA,dimABC:,dimAB:dimABC])
                        intsABCD.V[2][0][1][3] = numpy.array(Vraw[dimAB:dimABC,:dimA,dimA:dimAB,dimABC:])
                        intsABCD.V[3][0][1][2] = numpy.array(Vraw[dimABC:,:dimA,dimA:dimAB,dimAB:dimABC])
                        intsABCD.V[2][0][3][1] = numpy.array(Vraw[dimAB:dimABC,:dimA,dimABC:,dimA:dimAB])
                        intsABCD.V[3][0][2][1] = numpy.array(Vraw[dimABC:,:dimA,dimAB:dimABC,dimA:dimAB])
                        intsABCD.V[1][2][0][3] = numpy.array(Vraw[dimA:dimAB,dimAB:dimABC,:dimA,dimABC:])
                        intsABCD.V[1][3][0][2] = numpy.array(Vraw[dimA:dimAB,dimABC:,:dimA,dimAB:dimABC])
                        intsABCD.V[2][1][0][3] = numpy.array(Vraw[dimAB:dimABC,dimA:dimAB,:dimA,dimABC:])
                        intsABCD.V[3][1][0][2] = numpy.array(Vraw[dimABC:,dimA:dimAB,:dimA,dimAB:dimABC])
                        intsABCD.V[2][3][0][1] = numpy.array(Vraw[dimAB:dimABC,dimABC:,:dimA,dimA:dimAB])
                        intsABCD.V[3][2][0][1] = numpy.array(Vraw[dimABC:,dimAB:dimABC,:dimA,dimA:dimAB])
                        intsABCD.V[1][2][3][0] = numpy.array(Vraw[dimA:dimAB,dimAB:dimABC,dimABC:,:dimA])
                        intsABCD.V[1][3][2][0] = numpy.array(Vraw[dimA:dimAB,dimABC:,dimAB:dimABC,:dimA])
                        intsABCD.V[2][1][3][0] = numpy.array(Vraw[dimAB:dimABC,dimA:dimAB,dimABC:,:dimA])
                        intsABCD.V[3][1][2][0] = numpy.array(Vraw[dimABC:,dimA:dimAB,dimAB:dimABC,:dimA])
                        intsABCD.V[2][3][1][0] = numpy.array(Vraw[dimAB:dimABC,dimABC:,dimA:dimAB,:dimA])
                        intsABCD.V[3][2][1][0] = numpy.array(Vraw[dimABC:,dimAB:dimABC,dimA:dimAB,:dimA])
                        intsABCd += [intsABCD]
                intsABcd += [intsABCd]
            intsAbcd + [intsABcd]
        ints4 += [intsAbcd]
    return ints4
