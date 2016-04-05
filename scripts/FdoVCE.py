#!/usr/bin/python
"""
To run as a script:

    $ FdoVCE.py Aeff4 [Aeff5 Aeff6] outfile [he4fname he5fname he6fname]

Generates an interaction file based on He4, He5, and He6 output files.

If Aeff5 and Aeff6 are not specified, asssumes they are the same as Aeff4.
If filenames for he 4, 5, and 6 are not specified, assumes the convention
he[i]_[Aeff_i]/he[i]_A[Aeff_i].out is used, in which [i] is 4, 5, 6 and
[Aeff_i] is Aeff4, Aeff5, Aeff6.
"""

from sys import argv
from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException


def get_e0(aeff, fname_he4):
    if fname_he4 is not None:
        fname = fname_he4
    else:
        fname = 'he4_%d/he4_A%d.out' % (aeff, aeff)
    f = open(fname)
    line = f.readline()
    while 'State # 1   Energy' not in line:
        line = f.readline()
    f.close()
    ldat = line.split()
    e0 = float(ldat[5])
    return e0


def get_spe(aeff, e0, fname_he5):
    np1, np3 = 999., 999.
    if fname_he5 is not None:
        fname = fname_he5
    else:
        fname = 'he5_%d/he5_A%d.out' % (aeff, aeff)
    f = open(fname)
    for line in f:
        if 'State # ' not in line:
            continue
        ldat = line.split()
        e = float(ldat[5])
        j = int(2 * float(ldat[8]))
        if j == 1 and np1 == 999.:
            np1 = e - e0
        if j == 3 and np3 == 999.:
            np3 = e - e0
        if np1 != 999. and np3 != 999.:
            break
    f.close()
    return np1, np3


def print_header(aeff, e0, spe, st='%d'):
    header_lines = list()
    header_lines.append(
        '!  Effective SM interaction generated by OLS and VCE with Aeff = ' +
        st % (aeff,))
    header_lines.append('!  Zero body term: %10.6f' % (e0,))
    header_lines.append('!  Index  n  l  j tz')
    header_lines.append('!  1     0  1  1  1')
    header_lines.append('!  2     0  1  3  1')
    header_lines.append('! ')
    header_lines.append(
        '-999 ' + '  '.join(['%10.6f' % e for e in spe]) + '  4  6  0.000000')
    return '\n'.join(header_lines)


def get_tbme(aeff, e0, spe, fname_out, fname_he6, presc=None):
    write_lines = list()
    if presc is not None and (presc[0] != aeff or presc[1] != aeff):
        write_lines.append(print_header(str(presc), e0, spe, '%s'))
    else:
        write_lines.append(print_header(aeff, e0, spe))
    if fname_he6 is not None:
        fname = fname_he6
    else:
        fname = 'he6_%d/Heff_OLS.dat' % aeff
    f = open(fname)
    line = f.readline()
    dim = int(line.split()[0])
    kets = []
    for i in range(dim):
        ldat = f.readline().split()
        p, q = [int(dat) for dat in ldat[1:3]]
        j, t = [int(dat) for dat in ldat[9:11]]
        kets.append({'p': p, 'q': q, 'J': j, 'T': t})
    for i in range(dim):
        ldat = f.readline().split()
        for j in range(i, dim):
            if (kets[i]['J'], kets[i]['T']) != (kets[j]['J'], kets[j]['T']):
                continue
            v = float(ldat[j])
            if i == j:
                v -= e0 + spe[kets[i]['p'] - 1] + spe[kets[i]['q'] - 1]
            next_line = '%3d %3d %3d %3d  %3d %3d  %10.6f' % (
                kets[i]['p'], kets[i]['q'], kets[j]['p'], kets[j]['q'],
                kets[i]['J'], kets[i]['T'], v)
            write_lines.append(next_line)
    f.close()
    # write the file
    outfile = open(fname_out, 'w')
    outfile.write('\n'.join(write_lines))
    outfile.close()


def run(presc, fname_out, fname_he4, fname_he5, fname_he6):
    e0 = get_e0(presc[0], fname_he4)
    spe = get_spe(presc[1], e0, fname_he5)
    get_tbme(presc[2], e0, spe, fname_out, fname_he6, presc)

    # Do the inconsistent/universal way
    # Aeff = 6
    # E0 = GetE0(Aeff-2)
    # SPE = GetSPE(Aeff-1,E0)
    # GetTBME(Aeff,E0,SPE)


if __name__ == "__main__":
    if len(argv) == 3:
        a_prescription = (int(argv[1]),)*3
        out_fname = argv[2]
        he4_fname, he5_fname, he6_fname = (None,)*3
    if len(argv) == 6:
        a_prescription = (int(argv[1]),)*3
        out_fname, he4_fname, he5_fname, he6_fname = argv[2:6]
    elif len(argv) == 8:
        a_prescription = tuple([int(x) for x in argv[1:4]])
        out_fname = argv[4]
        he4_fname = argv[5]
        he5_fname = argv[6]
        he6_fname = argv[7]
    else:
        raise InvalidNumberOfArgumentsException(
            '\nFdoVCE.py called with %d arguments. ' % (len(argv)-1,) +
            'Please call with 2, 5, or 7 arguments.\n'
        )
    run(presc=a_prescription, fname_out=out_fname,
        fname_he4=he4_fname, fname_he5=he5_fname, fname_he6=he6_fname)
