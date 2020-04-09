#!/usr/bin/env python

from functions import *
import copy

Order = 3
diags, syms, bases = get_diags(Order)

dnum = 0
with open("DiagPolar{0}.new.txt".format(Order), 'a') as f:
    f.write("# Order {0} \n".format(Order))
for i in range(len(diags)):
#    if i != 7:
#        continue
    permutation = list(diags[i])
    print list(diags[i])
    order = int(len(permutation) / 2 + 1)
    sym = 1.0 / syms[i]
    loop_bases = bases[i]
    cancel_loops = find_loops_Hugen(permutation)
    Hugenholtz_diags = add_vertex(permutation, loop_bases, sym)
    # Remove the dressed Green and Hartree-Fock diagrams.
    for key, value in Hugenholtz_diags.items():
        if has_dressed_green(key, value[0]) or HasFock(key) or HasHartree(key):
            del Hugenholtz_diags[key]
        # elif has_dressed_green(key, value[0]):
        #     p = Diag(list(key), value[1], value[0])
        #     Diag.name = Diag.name + 1
        #     write(p)
        #     del Hugenholtz_diags[key]

    # Pick the diagrams which can maximize the cancellation.
    tem = copy.deepcopy(Hugenholtz_diags)
    topology_class = Group(tem)
    topology_num = len(topology_class)
    if topology_num == 0:
        continue
    num_list = [len(k) for k in topology_class]
    sym_list = [len(k)*1.0*sym for k in topology_class]
    required_num = [int(1.0*i/min(num_list)) for i in num_list]
    fin_diags = [[] for i in range(topology_num)]
    cancel_loops_list = []
    cancel_loops_poslist = []
    for num in cancel_loops.keys():
        cancel_loops_list += cancel_loops[num]
    loops_used = [cancel_loops_list[0]]

    for loop_num in range(1, len(cancel_loops_list)):
        # print 1, cancel_loops_list[loop_num], loops_used
        for t0 in [permutation[p] for p, v in enumerate(cancel_loops_list[loop_num]) if v == 1]:  # +[1]:
            for insert_1 in loops_used:
                for t1 in [permutation[p] for p,v in enumerate(insert_1) if v==1]:  # +[0]:
                    if (t0==t1): # or (t0==1) or (t1==0):
                        continue
                    tpgy, per = find_t0t1_diag([t0,t1], topology_class)
                    if (tpgy > -1) and (len(fin_diags[tpgy])<required_num[tpgy]):
                        if per not in fin_diags[tpgy]:
                            fin_diags[tpgy].append(per)
                    if all([len(fin_diags[i]) >= required_num[i] for i in range(len(required_num))]):
                        break
        # print 2, loops_used, cancel_loops_list[loop_num]
        for t1 in [permutation[p] for p, v in enumerate(cancel_loops_list[loop_num]) if v == 1]:  # +[0]:
            for insert_0 in loops_used:
                for t0 in [permutation[p] for p, v in enumerate(insert_0) if v == 1]:  # +[1]:
                    if t0 == t1:
                        continue
                    tpgy, per = find_t0t1_diag([t0,t1], topology_class)
                    if (tpgy > -1) and (len(fin_diags[tpgy])<required_num[tpgy]):
                        if per not in fin_diags[tpgy]:
                            fin_diags[tpgy].append(per)
                    if all([len(fin_diags[i]) >= required_num[i] for i in range(len(required_num))]):
                        break
        loops_used = loops_used + [cancel_loops_list[loop_num]]
    for loop_num in range(len(cancel_loops_list)):
        # print 3, cancel_loops_list[loop_num], cancel_loops_list[loop_num]
        for t0 in [permutation[p] for p, v in enumerate(cancel_loops_list[loop_num]) if v == 1]:  # +[1]:
            for t1 in [permutation[p] for p, v in enumerate(cancel_loops_list[loop_num]) if v == 1]:  # +[0]:
                if t0 == t1:
                    continue
                tpgy, per = find_t0t1_diag([t0, t1], topology_class)
                if (tpgy > -1) and (len(fin_diags[tpgy]) < required_num[tpgy]):
                    if per not in fin_diags[tpgy]:
                        fin_diags[tpgy].append(per)
                if all([len(fin_diags[i]) >= required_num[i] for i in range(len(required_num))]):
                    break
    sym_list = [1.0*sym_list[k]/len(fin_diags[k]) for k in range(len(fin_diags))]
    for tpgy in range(len(fin_diags)):
        for per in fin_diags[tpgy]:
            sym_factor = sym_list[tpgy]
            bases1 = Hugenholtz_diags[per][0]
            p = Diag(list(per), sym_factor, bases1)
            Diag.name = Diag.name + 1
            dnum = dnum + 1
            write(p)

print dnum

