#!/usr/bin/python
"""FdoVCE.py
Script of Ragnar Stroberg, edited for compatibility with my code.

To run as a script:

    $ FdoVCE.py Aeff4 Aeff5 Aeff6 A4 A5 A6 outfile
    [he4fname he5fname he6fname [nshell]]

Generates an interaction file based on He4, He5, and He6 output files.
"""

from sys import argv
from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException


class GroundStateEnergyNotFoundException(Exception):
    pass


class SingleParticleEnergyNotFoundException(Exception):
    pass


def get_j2_range(nshell):
    """Given the shell and whether the reversed convention is used, returns
    the ordered list of 2*j values for which single particle energies should
    be retrieved
    :param nshell: 0=s, 1=p, 2=sd, ...
    :return j2_range, idx_range
    j2_range is the ordered list of 2*j values
    idx_range is the ordered list of the indices of those 2*j values with
    respect to increasing j.
    For example, for
        j2_range  = [7/2, 1/2, 3/2, 5/2]
        idx_range = [  4,   1,   2,   3]
    """
    j2_i_range = [(1 + 2*xx, xx + 1) for xx in range(nshell + 1)]
    if nshell == 2:
        j2_i_range = list(reversed(j2_i_range))
    elif nshell != 1:
        print (
            'The correct convention for ordering SPEs is not known for '
            'nshell %d. They will be ordered by increasing j' % nshell
        )
    j2_range = [tup[0] for tup in j2_i_range]
    idx_range = [tup[1] for tup in j2_i_range]
    return j2_range, idx_range


def get_e0(fpath):
    """Given the file path to the first NCSD outfile (helium 4 for Nshell=1),
    returns the first energy with J=0
    :param fpath: file path to the first NCSD outfile (for Nshell=1, this
    would be helium 4)
    :return: first energy whose angular momentum is 0.0
    :raises GroundStateEnergyNotFoundException: raised when a state with J=0
    is not found in the given file
    """
    f = open(fpath)
    for line in f:
        if 'State #' in line:
            ldat = line.split()
            j = float(ldat[8])
            if j == 0.0:
                return float(ldat[5])
    else:
        raise GroundStateEnergyNotFoundException(
            '\nA ground state could not be retrieved from %s' % fpath)


def get_spe(fpath, e0, j2_range):
    """Given the path to the second NCSD outfile, returns a list of
    single particle energies ordered based on the ordering of j2_range
    :param fpath: path to the second NCSD outfile (helium 5 for Nshell=1)
    :param e0: zero body term (retrieved by get_e0)
    :param j2_range: ordered list of 2*j values for which SPE's exist for
    the shell. This should be listed in the order SPE's are to be written
    in the interaction file
    """
    spe_list = [999.] * len(j2_range)
    f = open(fpath)
    for line in f:
        if 'State # ' not in line:
            continue
        ldat = line.split()
        e = float(ldat[5])
        j = int(2 * float(ldat[8]))
        if j in j2_range:
            ii = j2_range.index(j)
            if spe_list[ii] == 999.:
                spe_list[ii] = e - e0
        if 999. not in spe_list:
            break
    else:
        raise SingleParticleEnergyNotFoundException(
            '\nOne or more SPE(s) could not be retrieved from %s' % fpath)
    f.close()
    return spe_list


def get_header_string(aeff_str, e0, spe, j2_range, a_values):
    """Returns the header to the interaction file, along with the zero body
    term and single particle energies
    :param aeff_str: Aeff for header line
    :param e0: zero body term
    :param spe: list of single particle energies ordered by increasing j
    :param j2_range: list of 2*j values for single particle energies, in the
    order they are to be printed
    :param a_values: first 3 exact A values corresponding to the A-prescription
    :return: header lines
    """
    header_lines = list()
    header_lines.append(
        '!  Effective SM interaction generated by OLS and VCE with Aeff = ' +
        '%s' % str(aeff_str))
    header_lines.append('!  Zero body term: %10.6f' % (e0,))
    header_lines.append('!  Index  n  l  j tz')
    for j, idx in zip(j2_range, range(len(spe))):
        header_lines.append('!  %d     %d  %d  %d  %d' % (idx+1, 0, 1, j, 1))
    header_lines.append('! ')
    spe_line = '-999 ' + '  '.join(['%10.6f' % e for e in spe])
    spe_line += '  %d  %d  0.000000' % (a_values[0], a_values[2])
    header_lines.append(spe_line)
    return '\n'.join(header_lines)


def get_tbme(aeff, e0, spe, j2_range, idx_range, a_values,
             fpath_write_int, fpath_heff, presc=None):
    """Writes interaction file based on Valence Cluster Expansion
    :param aeff: Aeff used for 3rd NCSD (helium6 if Nshell=1)
    :param e0: core energy
    :param spe: list of single particle energies (in order of increasing j)
    :param j2_range: ordered list of 2*j values for which SPE's exist for
    the shell. This should be listed in the order SPE's are to be written
    in the interaction file.
    :param idx_range: list of indices numbered according to increasing j and 
    ordered according to the same convention as j2_range.
    For example if the ordering of j's was [j=7/2, j=1/2, j=3/2, j=5/2],
    idx_range would be [4, 1, 2, 3].
    :param a_values: first 3 exact A values corresponding to the A-prescription
    :param fpath_write_int: file path for NuShellX interaction file to write
    :param fpath_heff: file path of Heff_OLS matrix
    :param presc: 3-tuple representing the A-prescription. If None, assumed to
    be (aeff, aeff, aeff)
    """
    # get the top header lines to write to *.int file
    write_lines = list()
    if presc is not None and (presc[0] != aeff or presc[1] != aeff):
        aeff_str = str(presc)
    else:
        aeff_str = str(aeff)
    write_lines.append(get_header_string(
        aeff_str=aeff_str, e0=e0, spe=spe, j2_range=j2_range, a_values=a_values
    ))
    # open Heff_OLS file for reading
    f = open(fpath_heff)
    line = f.readline()
    dim = int(line.split()[0])
    kets = []
    for i in range(dim):
        ldat = f.readline().split()
        p, q = [int(dat) for dat in ldat[1:3]]
        p = idx_range.index(p) + 1
        q = idx_range.index(q) + 1
        j, t = [int(dat) for dat in ldat[9:11]]
        kets.append({'p': p, 'q': q, 'J': j, 'T': t})
    # get TBMEs
    tbme_list = list()
    for i in range(dim):
        ldat = f.readline().split()
        for j in range(i, dim):
            if (kets[i]['J'], kets[i]['T']) != (kets[j]['J'], kets[j]['T']):
                continue
            v = float(ldat[j])
            if i == j:
                v -= (e0 + spe[kets[i]['p'] - 1] + spe[kets[i]['q'] - 1])
            tbme_list.append((
                kets[i]['p'], kets[i]['q'], kets[j]['p'], kets[j]['q'],
                kets[i]['J'], kets[i]['T'], v))
    f.close()
    # write TBME's in general convention:
    # a <= b, a <= c <= d
    next_tbme_list = list()
    for a, b, c, d, j, t, v in tbme_list:
        if a > b:
            a, b = b, a
        if c > d:
            c, d = d, c
        if a > c:
            a, b, c, d = c, d, a, b
        next_tbme_list.append((a, b, c, d, j, t, v))
    # sort TBME's by j, t, a, b, c, d
    next_tbme_list = sorted(next_tbme_list, key=lambda e: (e[4], e[5], e[:4]))
    # make lines
    for tbme in next_tbme_list:
        next_line = '%3d %3d %3d %3d  %3d %3d  %10.6f' % tbme
        write_lines.append(next_line)
    # write the file
    outfile = open(fpath_write_int, 'w')
    outfile.write('\n'.join(write_lines))
    outfile.close()


def run(presc, fpath_write_int, fpath_he4, fpath_he5, fpath_heff_ols, nshell,
        a_values):
    """Do the full valence cluster expansion based on the given A-prescription
    and NCSD file paths
    :param presc: 3-tuple representing the A-prescription
    :param fpath_write_int: path to the NuShellX interaction file to write
    :param fpath_he4: path to the 1st NCSD out file (helium 4 if Nshell=1)
    :param fpath_he5: path to the 2nd NCSD out file (helium 5 if Nshell=1)
    :param fpath_heff_ols: path to the Heff file generated by TRDENS
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd,...)
    :param a_values: first 3 exact A values
    """
    e0 = get_e0(fpath=fpath_he4)
    j2_range, idx_range = get_j2_range(nshell=nshell)
    spe = get_spe(fpath=fpath_he5, e0=e0, j2_range=j2_range)
    get_tbme(
        aeff=presc[2], e0=e0, spe=spe, presc=presc,
        j2_range=j2_range, idx_range=idx_range,
        fpath_write_int=fpath_write_int, fpath_heff=fpath_heff_ols,
        a_values=a_values,
    )

    # Do the inconsistent/universal way
    # Aeff = 6
    # E0 = GetE0(Aeff-2)
    # SPE = GetSPE(Aeff-1,E0)
    # GetTBME(Aeff,E0,SPE)


if __name__ == "__main__":
    nshell0 = 1
    if len(argv) == 8:
        a_prescription0 = tuple([int(x) for x in argv[1:4]])
        a_values0 = tuple([int(x) for x in argv[4:7]])
        out_fname = argv[7]
    elif len(argv) == 11:
        a_prescription0 = tuple([int(x) for x in argv[1:4]])
        a_values0 = tuple([int(x) for x in argv[4:7]])
        out_fname = argv[7]
        he4_fname = argv[8]
        he5_fname = argv[9]
        he6_fname = argv[10]
    elif len(argv) == 12:
        a_prescription0 = tuple([int(x) for x in argv[1:4]])
        a_values0 = tuple([int(x) for x in argv[4:7]])
        out_fname = argv[7]
        he4_fname = argv[8]
        he5_fname = argv[9]
        he6_fname = argv[10]
        nshell0 = argv[11]
    else:
        raise InvalidNumberOfArgumentsException(
            '\nFdoVCE.py called with %d arguments. ' % (len(argv)-1,) +
            'Please call with 7, 10, or 11 arguments.\n'
        )
    # TODO: make this script behavior better
    run(presc=a_prescription0, a_values=a_values0, fpath_write_int=out_fname,
        fpath_he4=he4_fname, fpath_he5=he5_fname, fpath_heff_ols=he6_fname,
        nshell=nshell0)
