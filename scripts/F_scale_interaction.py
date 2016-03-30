#!/usr/bin/python
"""F_scale_interaction.py

Reads the interaction TBME file from the given src file path, and writes it to
dst, changing the float values in accordance with scalefn.
"""

from os import path
from sys import argv
from ncsm_vce_calc import InvalidNumberOfArgumentsException

SCALE_FACTOR = 1.0
_LINE_FMT = ' %2d' * 6 + ' %11.6f' * 6 + '\n'


def ecf_off_diag_outside_valence(a_max):
    """Returns an elt_cond_fn that returns true if matrix element is off
    diagonal and outside the valence range
    :param a_max: maximum value in the valence space
    """
    def elt_cond_fn(elt):
        a, b, c, d, j, t, = elt
        if (a > a_max or b > a_max or c > a_max or d > a_max and
                set([a, b]) != set([c, d])):
            return True
        else:
            return False
    return elt_cond_fn


def ecf_isospin(t):
    """Returns an elt_cond_fn that returns true (do scaling) if isospin is
    the given value
    :param t: isospin
    :return: elt_cond_fn
    """
    def elt_cond_fn(elt):
        return elt[5] == t
    return elt_cond_fn


def get_scaleif_fn(elt_cond_fn, rest_cols, scalefactor):
    """Returns a scale function that does the scaling only if
    elt_cond_fn is satisfied
    :param elt_cond_fn: function from the tuple (a, b, c, d, j, t) to a
    boolean value
    :param rest_cols: columns of float values (relative to the first) to be
    scaled
    :param scalefactor: float value by which to scale rest_cols
    :return: function that takes elt and rest as arguments and returns a
    tuple of scaled "rest" values, where
    elt is the tuple (a, b, c, d, j, t) and
    rest is the tuple of associated float values,
    """
    def scaleif(elt, rest):
        if elt_cond_fn(elt):
            for col in rest_cols:
                rest[col] = scalefactor * rest[col]
        return rest
    return scaleif


def run(fpath_src, fpath_dst, scalefn0=lambda a, b: b, force=False,
        _line_fmt=_LINE_FMT):
    """Read and interaction file from fpath_src, scale it according to
    scalefn0, and write the resultant file into fpath_dst
    :param fpath_src: path to read
    :param fpath_dst: path to write
    :param scalefn0: function from two tuples to a tuple of scaled values.
    The first tuple is (a, b, c, d, j, t)
    The second tuple is the remaining columns
    The returned tuple is a scaled version of the second tuple
    :param force: if true, will overwrite an existing file at the destination
    path
    :param _line_fmt: format string for a line
    """
    if force or not path.exists(fpath_dst):
        fin = open(fpath_src, 'r')
        fout = open(fpath_dst, 'w')
        first_line = fin.readline()
        fout.write(first_line)
        for line in fin:
            ldat = line.split()
            elt = [int(x) for x in ldat[:6]]  # a, b, c, d, j, t
            rest = [float(x) for x in ldat[6:]]
            next_rest = scalefn0(elt, rest)
            next_line = _line_fmt % (tuple(elt) + tuple(next_rest))
            fout.write(next_line)
        fout.close()
        fin.close()
