#!/usr/bin/env python

import copy

def normalize(l):
    # normalize the cycle, pick up the cycle which the first element is maximum.
    return max([l[i:]+l[:i] for i in range(len(l))])


def even(x):
    if x < 2:
        return x
    if x % 2 == 0:
        return x
    else:
        return x - 1


def mirror(index):
    if index%2 == 0:
        return index + 1
    else:
        return index - 1


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
    print loops
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

permutation = [6,5,2,7,4,3]
print find_loops_Hugen(permutation)
