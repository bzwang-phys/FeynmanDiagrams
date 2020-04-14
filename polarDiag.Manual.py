#!/usr/bin/env python3

import sys, os
import numpy as np
from permutation import Permutation


def mirror(index):
    if index%2 == 0:
        return index + 1
    else:
        return index - 1


def add_vertex(Diag, ZMomentum, Sym):
    PolarDict = dict()
    SHIFT = 2
    order = int(len(Diag)/2+1)
    for i in range(2, len(Diag)+2):
        #Initialization
        #d[i]<== 1 <== 0 <== i
        d = [0,1]+list(Diag)
        d[1] = d[i]
        d[0] = 1
        d[i] = 0

        Momentum = np.zeros([order+1, 2*order], dtype=int)
        Momentum[1:, 2:] = ZMomentum
        Momentum[1:, 0] = ZMomentum[:, i-SHIFT]
        Momentum[1:, 1] = ZMomentum[:, i-SHIFT]
        Momentum[0, 0] = 1

        # if not CheckConservation(d, Momentum, IsPolarization=True):
        #     print("Momentum does not conserve or rank is wrong!")
        #     sys.exit(0)

        # print "Start with: ", d
        PolarDict[tuple(d)] = [Momentum, Sym]
        ToVisit = [d[1], mirror(d[1])]
        StartPermu = [tuple(d), tuple(d)]
        StartMom = [Momentum, Momentum]
        Visited = [0]
        while len(ToVisit) != 0:
            Index = ToVisit.pop()
            Permutation = list(StartPermu.pop())

            Mom = np.copy(StartMom.pop())
            if Index in Visited:
                continue

            if Permutation[1]!=Index and Permutation[1]!=mirror(Index):
                print("wrong!", Permutation, Index)
                sys.exit()
            Target = Permutation[Index]
            NextVertex = Permutation[1]
            PrevVertex = Permutation.index(1)
            Permutation[1] = Target
            Permutation[PrevVertex] = NextVertex
            Permutation[Index] = 1

            deltaMom = np.copy(Mom[:,PrevVertex]-Mom[:,1])
            Mom[:,1] = Mom[:,Index]
            Mom[:,Index] += deltaMom

            # if not CheckConservation(Permutation, Mom, IsPolarization=True):
            #     print "Momentum does not conserve or rank is wrong!"
            #     sys.exit(0)

            PolarDict[tuple(Permutation)] = [Mom, Sym]
            Visited.append(Index)

            if Target not in Visited:
                ToVisit.append(Target)
                ToVisit.append(mirror(Target))
                StartPermu.append(tuple(Permutation))
                StartPermu.append(tuple(Permutation))
                StartMom.append(Mom)
                StartMom.append(Mom)
        # print len(Visited)
    return PolarDict


def ver_loopbases(per, bases):
    points = len(per)
    order = int(points/2)
    ver_out_index = [i for i in range(points)]
    ver_in_index = [per.index(i) for i in range(points)]
    ver_bases = np.zeros((order+1, 2*(order-1)), dtype=int)
    for i in range(2, points, 2):
        ver_bases[:, i-2] = - bases[:, ver_out_index[i]] + bases[:, ver_in_index[i]]
        ver_bases[:, i-1] = - bases[:, ver_out_index[i+1]] + bases[:, ver_in_index[i]]
    return ver_bases


class Diag():
    name = 0

    def __init__(self, per, sym_factor, bases):
        self.sym_factor = sym_factor
        self.permutation = per
        self.loop_bases = bases
        self.ver_bases = ver_loopbases(self.permutation, self.loop_bases)


def write(P):
    # --------------- symmetry factor ----------
    if order == 4:
        a1, a2, a3, a4, a5, a6, a7, a8 = P.permutation
        pm = Permutation(a1+1, a2+1, a3+1, a4+1, a5+1, a6+1, a7+1, a8+1)
    elif order == 3:
        a1, a2, a3, a4, a5, a6 = P.permutation
        pm = Permutation(a1+1, a2+1, a3+1, a4+1, a5+1, a6+1)
    sign = pm.sign
    P.sym_factor = -1 * sign * abs(P.sym_factor)
    # ---------------------------------------------
    with open("DiagPolar{0}.new.txt".format(order), 'a') as f:
        f.write("# Topology   [" + str(P.name) + "]\n")
        for i in P.permutation:
            f.write('{0:3d}'.format(i))
        f.write('\n')
        f.write("# Symmetry Factor \n {} \n".format(P.sym_factor))
        f.write("# Loop Bases" + '\n')
        len_i, len_j = P.loop_bases.shape
        for i in range(len_i):
            for j in range(len_j):
                f.write('{0:3d}'.format(P.loop_bases[i, j]))
            f.write('\n')
        f.write("# Ver Loop Bases" + '\n')
        len_i, len_j = P.ver_bases.shape
        for i in range(len_i):
            for j in range(len_j):
                f.write('{0:3d}'.format(P.ver_bases[i, j]))
            f.write('\n')
        f.write('\n')


permutation = [4, 5, 6, 3, 2, 7]
loop_bases = np.array([[1,0,0,1,0,0], [0,1,0,1,0,0], [0,1,1,0,1,0],[0,0,0,0,0,1]], dtype=int)
# permutation = [4, 5, 2, 3]
# loop_bases = np.array([[1,0,1,0], [0,1,1,0], [0,1,0,1]], dtype=int)
sym = -1.0
all_diag = add_vertex(permutation, loop_bases, sym)
order = int(len(permutation) / 2 + 1)

with open("DiagPolar{0}.new.txt".format(order), 'a') as f:
    f.write("# Order {0} \n".format(order))
for key, value in all_diag.items():
    if (list(key)[0] == 1) or (list(key)[1] == 0):
        continue
    if (list(key)[0] == 7) and ((list(key)[1] == 3) or (list(key)[1] == 3)):
        p = Diag(list(key), value[1], value[0])
        Diag.name = Diag.name + 1
        write(p)


