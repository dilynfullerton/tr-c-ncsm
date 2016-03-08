#!/usr/bin/python
"""
$ python GetSmallerInteraction.py <fname> <n1> <n2> <outfile_name>

Reads the TBME file and removes elements beyond the
cutoff given by n1. 
NOTE: These are printed out, so the
result should be piped to the desired file. Then,
the final line needs to be moved to the top.
"""

from sys import argv
from os import path


def get_n(indx):
    n = 0
    while (n+1)*(n+2)/2 < indx: 
        n += 1
    return n


def run(fname, n1, n2, f_out_name=None):
    if f_out_name is not None and path.exists(f_out_name):
        # test whether the operation has already been done
        f = open(f_out_name, 'r')
        n1_bef, n2_bef = [int(x) for x in f.readline().split()[1:3]] 
        if n1_bef == n1 and n2_bef == n2:
            return 0
    f = open(fname)
    
    line = f.readline()
    ldat = line.split()
    hw, a_tbme = [float(x) for x in ldat[3:]]

    index_max = (n1+1)*(n1+2)/2
    
    next_lines = list()
    ntbme = 0
    for line in f:
        ldat = line.split()
        a, b, c, d, j, t = [int(x) for x in ldat[:6]]
        if a > index_max or b > index_max or c > index_max or d > index_max:
            continue
        if get_n(a)+get_n(b) > n2 or get_n(c)+get_n(d) > n2:
            continue
        next_lines.append(line.rstrip('\n'))
        ntbme += 1
    
    s = '%d  %d  %d  %.4f  %.4f' % (ntbme, n1, n2, hw, a_tbme)
    next_lines.insert(0, s)
    f.close()
    
    if f_out_name is not None:
        f_out = open(f_out_name, 'w')
        for line in next_lines:
            f_out.write(line)
            f_out.write('\n')
        f_out.close()
    else:
        for line in next_lines:
            print line
    return 1


# if len(argv) >= 4:
#     fname0 = argv[1]
#     n1_0, n2_0 = [int(x) for x in argv[2:4]]
#     if len(argv) >= 5:
#         f_out_name0 = argv[4]
#     else:
#         f_out_name0 = None
#     run(fname=fname0, n1=n1_0, n2=n2_0, f_out_name=f_out_name0)
