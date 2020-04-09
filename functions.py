#!/usr/bin/env python
import numpy as np
import unionfind
import random
from nullspace import rank, nullspace
from numpy.linalg import matrix_rank
import copy
import sys,os
from permutation import Permutation


random.seed(2)

def ReNameDiag(DiagListInput, Initial=1):
    DiagList=list(DiagListInput)
    for d in DiagList:
        for i in range(len(d)/2):
            d[i*2]=d[i*2]*2-2+Initial*2
            d[i*2+1]=d[i*2+1]*2+1-2+Initial*2
    return DiagList


def swap(array, i, j):
    array = list(array)
    array[i], array[j] = array[j], array[i]
    return tuple(array)

def IsConnected(permutation, reference, InteractionPairs):
    diagram=set(InteractionPairs)
    for i in range(len(permutation)):
        diagram.add((reference[i], permutation[i]))

    n_node = len(InteractionPairs)*2
    diagram_union = unionfind.UnionFind(n_node)

    for edge in diagram:
        if edge[0]!=edge[1] and not diagram_union.is_connected(edge[0], edge[1]):
            diagram_union.union(edge[0], edge[1])
    return diagram_union.get_n_circles() == 1


def GetInteractionPairs(Order):
    return [(2*i,2*i+1) for i in range(Order)]


def GetReference(Order):
    return range(2*Order)

def HasTadpole(permutation, reference):
    for i in range(len(permutation)):
        if reference[i]==permutation[i]:
            return True
    return False


def HasFock(permutation):
    for i in range(len(permutation)):
        # end=reference[i]
        end=permutation[i]
        if i==0 or i==1:
            continue
        if abs(i-end)==1 and min(i, end)%2==0:
            return True
    return False


def HasHartree(permutation):
    for i in range(len(permutation)):
        if i == permutation[i]:
            return True
    return False


def swap(array, i, j):
    array = list(array)
    array[i], array[j] = array[j], array[i]
    return tuple(array)

def swap_interaction(permutation, m, n, k, l):
    permutation = list(permutation)
    mp,np,kp,lp=(permutation.index(e) for e in (m,n,k,l))
    permutation[mp]=k
    permutation[kp]=m
    permutation[np]=l
    permutation[lp]=n
    permutation[m],permutation[k]=permutation[k],permutation[m]
    permutation[n],permutation[l]=permutation[l],permutation[n]
    return tuple(permutation)

def swap_LR(permutation, i, j):
    # print permutation, i, j
    permutation = list(permutation)
    ip,jp = permutation.index(i),permutation.index(j)
    permutation[ip] = j
    permutation[jp] = i
    permutation[i],permutation[j] = permutation[j],permutation[i]
    return tuple(permutation)

def swap_LR_Hugen(permutation, i, j):
    permutation = list(permutation)
    permutation[i],permutation[j] = permutation[j],permutation[i]
    return swap_LR(permutation, i, j)
    # return tuple(permutation)


def check_Unique_Permutation(permutation, PermutationDict, TimeRotation):
    Order = len(permutation)/2
    Deformation = [permutation]

    if TimeRotation:
        for idx in range(1, Order):
            for i in range(len(Deformation)):
                for j in range(1, idx):
                    Deformation.append(swap_interaction(Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

    for idx in range(1,Order):
        for i in range(len(Deformation)):
            Deformation.append(swap_LR(Deformation[i], idx*2, idx*2+1))

    for idx in range(1,Order):
        for i in range(len(Deformation)):
            Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

    Deformation = set(Deformation)
    DeformationFinal = []
    for p in Deformation:
        if p in PermutationDict:
            # DeformationFinal+=list(PermutationDict[p])
            del PermutationDict[p]
            DeformationFinal.append(p)

    # print "remaining length of permutation dictionary:", len(PermutationDict)
    return list(DeformationFinal)


def get_Unique_Permutation(permutationList, TimeRotation=True):
    Order = len(permutationList[0])/2
    PermutationDict={}
    for p in permutationList:
        PermutationDict[tuple(p)]=None
    for per in permutationList:
        if not PermutationDict.has_key(tuple(per)):
            continue
        Deformation = [per]

        if TimeRotation:
            for idx in range(1, Order):
                for i in range(len(Deformation)):
                    for j in range(1, idx):
                        Deformation.append(swap_interaction(Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

        for idx in range(1,Order):
            for i in range(len(Deformation)):
                Deformation.append(swap_LR(Deformation[i], idx*2, idx*2+1))

        # for idx in range(1,Order):
            # for i in range(len(Deformation)):
                # Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

        Deformation = set(Deformation)
        for p in Deformation:
            if tuple(p)==tuple(per):
                continue
            if p in permutationList:
                del PermutationDict[p]

    print "remaining length of permutation dictionary:", len(PermutationDict)
    return PermutationDict.keys()


def Group(PermutationDict, TimeRotation=True):
    UnlabeledDiagramList = []
    FactorList = []
    # for permutation in PermutationList[0:1]:
    while len(PermutationDict)>0:
        # print "Remaining diagram {0}".format(len(PermutationDict))
        permutation = PermutationDict.keys()[0]
        Deformation = check_Unique_Permutation(permutation, PermutationDict, TimeRotation)
        # if len(Deformation)>0:
        UnlabeledDiagramList.append(Deformation)
    return UnlabeledDiagramList


def FindIndependentK(permutation, reference, InteractionPairs):
    # kList=[(random.randint(0, Nmax), random.randint(0,Nmax)) for i in range(len(InteractionPairs)+1)]
    N=len(InteractionPairs)
    Matrix=np.zeros((2*N,3*N))
    for i in range(2*N):
        interaction=int(i/2)+2*N
        sign=i%2
        Matrix[i,interaction]=-(-1)**sign
        Matrix[i, i]=-1
        Matrix[i, permutation.index(i)]=1
    # print Matrix
    vectors = nullspace(Matrix)
    # print len(vectors)
    # print vectors
    freedoms=vectors.shape[1]
    if freedoms!=N+1:
        print "Warning! Rank is wrong for {0} with \n{1}".format(permutation, vectors)
    return vectors


def AssignMomentums(permutation, reference, InteractionPairs):
    N=len(InteractionPairs)
    vectors=FindIndependentK(permutation, reference, InteractionPairs)
    freedoms=vectors.shape[1]
    karray=np.array([random.random() for _ in range(freedoms)])
    kVector=np.dot(vectors, karray)
    # kVector=vectors[:,0]
    return kVector[:2*N], kVector[2*N:]


def GenerateAllDiagram(UnlabeledDiagram, InteractionPairs, DoesCheck=True):
    AllDiagram=[]
    for d in UnlabeledDiagram:
        Deformation=[d]
        for idx in range(1,Order):
            for i in range(len(Deformation)):
                Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

        DeformationFinal=list(Deformation)
        if DoesCheck is True:
            for p in Deformation:
                kG, kW=AssignMomentums(p, Reference, InteractionPairs)

                Flag=True
                for i in range(len(kW)):
                    if Flag and abs(kW[i])<1e-12:
                        # print "k=0 on W {0}: {1}".format(p, kW[i])
                        DeformationFinal.remove(p)
                        Flag=False
                        break

                for j in range(1,len(kW)):
                    if Flag and abs(abs(kW[0])-abs(kW[j]))<1e-12:
                        DeformationFinal.remove(p)
                        Flag=False
                        break

                for i in range(0,len(kG)):
                    for j in range(i+1,len(kG)):
                        if Flag and abs(kG[i]-kG[j])<1e-12:
                            # print "Same k on G for {0}: {1} on {2}; {3} on {4}".format(p, kG[i],i,kG[j],j)
                            # print "Same k on W for {0}: {1}; 1, {2}".format(p, kG[i],kG[j])
                            DeformationFinal.remove(p)
                            Flag=False
                            # print "Flag",Flag
                            break
        AllDiagram+=DeformationFinal
    return AllDiagram

def GenerateAllFreeEnergyDiagram(UnlabeledDiagram, InteractionPairs):
    Order = int(len(InteractionPairs))
    # print "Order", Order
    # print "Diagram", UnlabeledDiagram
    d = UnlabeledDiagram
    Deformation = [d]
    for idx in range(Order):
        for i in range(len(Deformation)):
            Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

    for idx in range(Order):
        for i in range(len(Deformation)):
            for j in range(idx):
                Deformation.append(swap_interaction(Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

    for idx in range(Order):
        for i in range(len(Deformation)):
            Deformation.append(swap_LR(Deformation[i], idx*2, idx*2+1))
    DeformationFinal=[tuple(e) for e in Deformation]
    return DeformationFinal


def StartPoint(Order):
    StartPoint=range(Order*2)
    Momentum=np.zeros([Order+1, 2*Order], dtype=int)
    Momentum[0,0]=1
    Momentum[-1,-1]=1
    for i in range(1,Order):
        StartPoint[i*2-1], StartPoint[i*2]=StartPoint[i*2], StartPoint[i*2-1]
        Momentum[i, i*2-1]=1
        Momentum[i, i*2]=1
    FermiSign = -1  #n+1 loop  contributes (-1)^(n+1) and order n contributes (-1)^n
    return tuple(StartPoint), Momentum, FermiSign


def CheckConservation(permutation, MomentumBases, IsPolarization = False):
    Order = len(permutation)/2
    if matrix_rank(MomentumBases) != Order+1:
        print "rank is wrong with permutation {0}\n{1}".format(permutation, MomentumBases)
        return False
    Momentum = np.zeros(2*Order)
    for i in range(Order):
        Momentum += random.random()*MomentumBases[i, :]
    for i in range(Order):
        In1, In2 = 2*i, 2*i+1
        Out1 = permutation.index(2*i)
        Out2 = permutation.index(2*i+1)
        TotalMom = Momentum[In1] + Momentum[In2] - Momentum[Out1] - Momentum[Out2]
        if abs(TotalMom)>1e-10:
            print "Vertex {0} breaks the conservation laws. Bases: \n{1}".format(i, Momentum)
            print In1, In2, Out1, Out2
            print permutation
            print MomentumBases
            return False
    if IsPolarization:
        #the first loop basis has to be the external momentum
        Ext = np.zeros(Order+1, dtype=int)
        Ext[0] = 1
        if not np.all(MomentumBases[:,0]-MomentumBases[:,permutation.index(0)]-Ext==0):
            print "The first loop basis is not the external momentum"
            print permutation, MomentumBases
            sys.exit(0)
        if not np.all(MomentumBases[:,1]-MomentumBases[:,permutation.index(1)]+Ext==0):
            print "The first loop basis is not the external momentum"
            print permutation, MomentumBases
            sys.exit(0)
    return True


def GenerateMomentum(permutation, OldMomentum, i, j):
    if i/2==j/2:
        return None
    Order=len(permutation)/2
    if permutation[i]/2==permutation[j]/2:
        Momentum=np.copy(OldMomentum)
    else:
        Momentum=np.copy(OldMomentum)
        ni=i/2
        nj=j/2
        ip=4*ni+1-i
        jp=4*nj+1-j
        Momentum[:,[i,j]]=Momentum[:,[j,i]]

        if permutation[ip]==jp or permutation[ip]==j:
            Momentum[:,ip]+=Momentum[:,j]-Momentum[:,i]
            # print "Connect ip to jp", i, j, ip, jp, permutation[ip], permutation[jp]
        elif permutation[jp]==i or permutation[jp]==ip:
            Momentum[:,jp]+=Momentum[:,i]-Momentum[:,j]
            # print "Connect jp to ip", i, j, ip, jp, permutation[ip], permutation[jp]
        else:
            return -1
    if not CheckConservation(permutation, Momentum):
        print "Conservation or Rank Check fails."
        sys.exit(0)
        return None
    return Momentum

def GetAllPermutations(Order, DiagDict, DiagInvDict, DiagSymDict):
    """ 
    output:
        Diagrams: a dictionary contains a map from the original diagram to the diagram with optimized bases
        OptDiagrams: a dictionary contains a map from the optimized diagram to a list (original diagram, momentum bases for the optimized diagram, the symmetry factor for the optimized diagram)
    """
    reference=GetReference(Order)
    InteractionPairs=GetInteractionPairs(Order)

    permutation, Momentum, FermiSign=StartPoint(Order)
    PermuList = [permutation]
    MomList=[Momentum]
    SignList=[FermiSign]

    idx = 0
    while idx < 2*Order:
        # print "Index {0}".format(idx)
        for i in range(len(PermuList)):
            for j in range(idx):
                newpermutation=tuple(swap(PermuList[i], idx, j))
                # print "Old {0}, New {1} by switching {3},{4}\n OldBases: \n{2}".format(PermuList[i], newpermutation, MomList[i], idx, j)
                if IsConnected(newpermutation, reference, InteractionPairs):
                    Momentum=GenerateMomentum(newpermutation, MomList[i], idx, j) 
                    if Momentum is not None and Momentum is not -1:
                        # print "old :{0}, old_bases: \n{1}\n new:{2}, switch {4},{5},  new_bases: \n {3}".format(PermuList[i], MomList[i], newpermutation, Momentum, idx, j)
                        PermuList.append(newpermutation)
                        MomList.append(Momentum)
                        SignList.append(-SignList[i])
        idx += 1

    Diagrams={}
    for i in range(len(PermuList)):
        p=PermuList[i]
        Diag=tuple(DiagInvDict[p])
        SymFactor=SignList[i]*abs(DiagSymDict[Diag])
        Diagrams[Diag]=(p, MomList[i], SymFactor)

    OptDiagrams={}
    for k in Diagrams.keys():
        # print "Diagram {0}: {1} with SymFactor {3}\n {2}".format(k, Diagrams[k][0], Diagrams[k][1], Diagrams[k][2])
        p, Mom, Sym=Diagrams[k]
        OptDiagrams[p]=(k, Mom, Sym)

    # print "Total Diagram {0} vs {1}".format(len(Diagrams.keys()),len(DiagDict.keys()))

    for k in DiagDict.keys():
        if Diagrams.has_key(k) is False:
            print k

    return Diagrams, OptDiagrams

def FindAllLoops(permutation):
    order = len(permutation)/2
    Visited = set()
    path=[]
    for e in permutation:
        newloop=[]
        vertex=e
        while vertex not in Visited:
            newloop.append(vertex)
            Visited.add(vertex)
            vertex=permutation[vertex]
        if len(newloop)>0:
            path.append(newloop)
    if sum([len(l) for l in path])!=2*order:
        print "length of all loops should be 2*order"
        sys.exit(0)
    return path

def SymmetrizeLoops(OptDiagrams, DiagDict, DiagInvDict):
    order = len(OptDiagrams.keys()[0])/2
    NewOptDiagrams = dict(OptDiagrams)
    for diag in OptDiagrams.keys():
        if diag not in NewOptDiagrams.keys():
            continue
        Momentum=OptDiagrams[diag][1]
        AllDiag=[list(diag)]
        AllMom=[Momentum]
        Loops=FindAllLoops(diag)
        # print "original",diag
        # print "loops",Loops
        for l in Loops:
            for k in range(len(AllDiag)):
                rDiag=list(AllDiag[k])
                rMom=np.array(AllMom[k])
                for e in l:
                    Next=diag[e]
                    rDiag[Next]=e
                    rMom[:,Next]=-Momentum[:,e]
                # print "generate", rDiag
                AllDiag.append(rDiag)
                AllMom.append(rMom)
        EqualDict={}
        OptDict={}
        for d,m in zip(AllDiag, AllMom):
            if not CheckConservation(d, m):
                print "Conservation check fails!"
                sys.exit(0)
            if DiagInvDict.has_key(tuple(d)):
                if not EqualDict.has_key(DiagInvDict[tuple(d)]):
                    EqualDict[DiagInvDict[tuple(d)]]=(d, m)
        print "{0} equals to {1}".format(diag, [EqualDict[k][0] for k in EqualDict.keys()])

        for k in NewOptDiagrams.keys():
            if k is diag:
                continue
            if EqualDict.has_key(DiagInvDict[tuple(k)]):
                d, m=EqualDict[DiagInvDict[tuple(k)]]
                value=NewOptDiagrams[k]
                new_value=(value[0], m, value[-1])
                NewOptDiagrams[tuple(d)]=new_value
                del NewOptDiagrams[k]

    print "Old", OptDiagrams.keys()
    print "New", NewOptDiagrams.keys()
    return NewOptDiagrams

# def GetWLoopBases(diag, mom):

def SaveToFile(UniqueDiagrams, Name):
    if len(UniqueDiagrams)==0: 
        return
    diag, sym, mom, all_diag=UniqueDiagrams[0]
    order=len(diag)/2
    with open("./Diag{0}{1}.txt".format(Name, order), "a") as f:
        f.write("#Order {0}\n".format(order))
        for diag, sym, mom, all_diag in UniqueDiagrams:
            f.write("# Topology\n")
            for i in diag:
                f.write("{0:2d} ".format(i))
            f.write("\n")
            f.write("# SymmetryFactor\n{0}\n".format(sym))
            f.write("# Loop Bases\n")
            for i in range(order+1):
                for j in range(2*order):
                    f.write("{0:2d} ".format(mom[i,j]))
                f.write("\n")

            print "all diagram to print for ", diag
            # print all_diag

            f.write("#Ver Loop Bases\n")


            ###### Generate indepdent Mom ################
            WMom=[]
            LoopNum=[]
            for j in range(1,order):
                end=2*j
                start=diag.index(end)
                WMom.append(mom[:,start]-mom[:, end])
                WMom.append(mom[:,start]-mom[:, end+1])

            for i in range(order+1):
                for j in range(2*(order-1)):
                    # end=2*j
                    # start=diag.index(end)
                    # print start, end, mom[i, start]-mom[i, end]
                    f.write("{0:2d} ".format(WMom[j][i]))
                f.write("\n")

            f.write("\n")

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
    order = len(P.permutation)/2
    # --------------- symmetry factor ----------
    if order == 6:
        a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12 = P.permutation
        pm = Permutation(a1+1, a2+1, a3+1, a4+1, a5+1, a6+1, a7+1, a8+1, a9+1, a10+1, a11+1, a12+1)
    elif order == 5:
        a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = P.permutation
        pm = Permutation(a1+1, a2+1, a3+1, a4+1, a5+1, a6+1, a7+1, a8+1, a9+1, a10+1)
    elif order == 4:
        a1, a2, a3, a4, a5, a6, a7, a8 = P.permutation
        pm = Permutation(a1+1, a2+1, a3+1, a4+1, a5+1, a6+1, a7+1, a8+1)
    elif order == 3:
        a1, a2, a3, a4, a5, a6 = P.permutation
        pm = Permutation(a1+1, a2+1, a3+1, a4+1, a5+1, a6+1)
    elif order == 2:
        a1, a2, a3, a4 = P.permutation
        pm = Permutation(a1+1, a2+1, a3+1, a4+1)
    sign = pm.sign
    P.sym_factor = sign * abs(P.sym_factor)
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


def get_diags(o):  
    Order = o
    Order -= 1  
    # global Order
    Reference=GetReference(Order)
    InteractionPairs=GetInteractionPairs(Order)

    FreeEnergyDiagramDict={}
    FreeEnergyDiagramSymDict={}
    FreeEnergyDiagramInvDict={}
    with open("./Diagram/HugenDiag{0}.txt".format(Order)) as f:
        d=f.read()
        exec(d)
        TotalSym=0
        DiagList=ReNameDiag(Diag, Initial=0)

        # print "lnZ diagram List:", DiagList
        # print "lnZ diagram Symmetry factor:", Sym

        for d, s in zip(DiagList, Sym):
            AllDiagrams=GenerateAllFreeEnergyDiagram(d, InteractionPairs)
            # print "{0} with SymFactor: {1}, and Number: {2} with duplicate {3}".format(d, s, len(AllDiagrams), len(AllDiagrams)/len(set(AllDiagrams))) 
            # print "{0}".format(set(AllDiagrams)) 
            TotalSym+=float(len(AllDiagrams))/abs(s)
            FreeEnergyDiagramDict[tuple(d)]=AllDiagrams
            FreeEnergyDiagramSymDict[tuple(d)]=s
            for e in AllDiagrams:
                FreeEnergyDiagramInvDict[tuple(e)]=tuple(d)

        # print "Total Free energy diagrams: {0}, TotalSym: {1}".format(len(FreeEnergyDiagramDict), TotalSym)


        DiagDict, OptDiagDict=GetAllPermutations(Order, FreeEnergyDiagramDict, FreeEnergyDiagramInvDict, FreeEnergyDiagramSymDict)

    Order+=1 #get the order for polarization
    Reference=GetReference(Order)
    InteractionPairs=GetInteractionPairs(Order)
    TotalDiagNum=0
    TotalSym=0.0

    TempDiagList=OptDiagDict.keys()
    DiagList=[]
    SymList=[]
    MomList=[]
    for d in TempDiagList:
        DiagList.append(tuple([e+2 for e in d]))
        MomList.append(OptDiagDict[d][1])
        SymList.append(OptDiagDict[d][2])

    return DiagList, SymList, MomList


def has_dressed_green(permutation, bases):
    length = len(permutation)
    flag = False
    for i in range(length):
        for j in range(i+1, length):
            if (bases[:,i] == bases[:,j]).all():
                flag = True
    return flag


def normalize(l):
    # normalize the cycle, pick up the cycle which the first element is maximum.
    return max([l[i:]+l[:i] for i in range(len(l))])


def find_loops_Hugen(permutation):
    order = int(len(permutation)/2)
    # find all loops in a directed graph (before add_vertex).
    def is_loop(num):
        if (num in trace) or (mirror(num) in trace):
            if num in trace:
                index = trace.index(num)
            elif mirror(num) in trace:
                index = trace.index(mirror(num))
            loops.append(trace[index:len(trace)+1])
            return
        trace.append(num)
        for j in range(even(num), even(num)+2):
            is_loop(permutation[j-2])
        trace.pop()
        return
    trace = []
    loops = []
    is_loop(2)
    is_loop(3)
    # Remove the duplicated loops.
    loops_unique = set(map(tuple, map(normalize, loops)))
    loops = map(list, list(loops_unique))
    # Rebuild the loop with permutation_org
    loops_dict = {}
    loops_dict_bk = {}
    list2pair = lambda x:[[x[i],x[i+1]] if i!=len(x)-1 else [x[i],x[0]] for i in range(len(x))]
    for loop in loops:
        length = len(loop)
        loops_dict_bk.setdefault(length, []).append(loop)
        loop = list2pair(loop)
        loop_new = [[0] * (2 * order)]
        for line in loop:
            if (permutation[line[0]-2] in [line[1], mirror(line[1])]) and (permutation[mirror(line[0])-2] not in [line[1], mirror(line[1])]):
                for ind in range(len(loop_new)):
                    loop_new[ind][line[0] - 2] = 1
            elif (permutation[line[0]-2] not in [line[1], mirror(line[1])]) and (permutation[mirror(line[0])-2] in [line[1], mirror(line[1])]):
                for ind in range(len(loop_new)):
                    loop_new[ind][mirror(line[0])-2] = 1
            elif (permutation[line[0] - 2] in [line[1], mirror(line[1])]) and (permutation[mirror(line[0]) - 2] in [line[1], mirror(line[1])]):
                tem = copy.deepcopy(loop_new)
                for ind in range(len(loop_new)):
                    loop_new[ind][line[0] - 2] = 1
                    tem[ind][mirror(line[0])-2] = 1
                loop_new.extend(tem)
        loops_dict.setdefault(length, []).extend(loop_new)
        loops_dict = {k:map(list, list(set(map(tuple, v)))) for k,v in loops_dict.items()}
    return loops_dict


def swap_Hugen_org(permutation, i, j):
    # description: exchange i and j in the permutation_org
    permutation = list(permutation)
    ip, jp = permutation.index(i), permutation.index(j)
    permutation[ip] = j
    permutation[jp] = i
    return permutation


def swap_Hugen_per(permutation, i, j):
    # description: exchange i and j in the permutation
    permutation = list(permutation)
    permutation[i], permutation[j] = permutation[j], permutation[i]
    return permutation


def topology_equal_diags(permutation):
    # description: return the list of diagrams which is topologically equivalent to the "permutation".
    order = len(permutation)/2
    topology_diags = [permutation]

    # All possible exchanges between Hugenholtz nodes.
    for i in range(len(topology_diags)):
        for jp in range(1, order):
            for kp in range(1, jp):
                topology_diags.append(swap_interaction(topology_diags[i], jp*2, jp*2+1, kp*2, kp*2+1))

    # All possible exchanges inside every Hugenholtz nodes.
    for i in range(len(topology_diags)):
        for idx in range(1, order):
            topology_diags.append(swap_Hugen_org(topology_diags[i], idx*2, idx*2+1))
    for i in range(len(topology_diags)):
        for idx in range(1, order):
            topology_diags.append(swap_Hugen_per(topology_diags[i], idx*2, idx*2+1))
    topology_diags = map(list,list(set(map(tuple,topology_diags))))
    return topology_diags


def topology_classify(diags):
    # description: Topological classification of diags, return a dictionary {1:[...], 2:[...],....}.
    topology_class = {}
    diags = map(list, diags)
    for i in range(len(diags)):
        if len(topology_class) == 0:
            topology_class[1] = [diags[i]]
            continue
        # There already have diags[i] in topology_class.
        flag = True
        for key, value in topology_class.items():
            if is_topological_equal(diags[i], value[0]):
                topology_class[key].append(diags[i])
                flag = False
                break
        # It doesn't have diags[i] in all_topology.
        if flag:
            index = max(topology_class.keys()) + 1
            topology_class[index] = [diags[i]]
    return topology_class


def even(x):
    if x < 2:
        return x
    if x % 2 == 0:
        return x
    else:
        return x - 1


def is_topological_equal(diag1, diag2):
    order = int(len(diag1) / 2)
    diag1 = [even(diag1[i]) for i in range(len(diag1))]  # permutation in even-representation
    diag2 = [even(diag2[i]) for i in range(len(diag2))]
    position1 = [[even(j) for j,k in enumerate(diag1) if k==find] for find in range(2,2*order,2)] # permutation_org in even-representation
    position2 = [[even(j) for j,k in enumerate(diag2) if k==find] for find in range(2,2*order,2)]

    position1.sort(key=lambda x:(x[0],x[1]))    # normalize the position_list,
    position2.sort(key=lambda x:(x[0],x[1]))    # use the position list to replace the permutation.
    res = (np.array(position1) == np.array(position2)).all()
    return res


def find_t0t1_diag(to, topo_class):
    for topo_index in range(len(topo_class)):
        for diag in topo_class[topo_index]:
            if diag[:2] == tuple(to):
                return topo_index, diag
    return -1, -1




