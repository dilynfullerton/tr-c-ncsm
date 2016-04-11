#!/usr/bin/python
"""F_scale_interaction.py

Reads the interaction TBME file from the given src file path, and writes it to
dst, changing the float values in accordance with scalefn.
"""

from __future__ import division
from os import path

SCALE_FACTOR = 1.0
_LINE_FMT = ' %3d' * 6 + ' %11.6f' * 6 + '\n'


def ecf_off_diag_outside_valence(nshell):
    """Returns an elt_cond_fn that returns true if matrix element is off
    diagonal and couples the valence + core space to the outside space
    or the core space to the valence space
    Assumes that the interaction file's format satisfies:
        (1) a <= b
        (2) a <= c
        (3) c <= d
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    """
    idx_max_core = nshell * (nshell + 1) / 2
    idx_max_valn = (nshell + 1) * (nshell + 2) / 2

    def elt_cond_fn(elt):
        a, b, c, d, j, t, = elt
        return ((b <= idx_max_valn < d) or (d <= idx_max_valn < b) or
                (b <= idx_max_core < d) or (d <= idx_max_core < b))
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


def get_scaleif_fn(elt_cond_fn, rest_cols_to_scale, scalefactor):
    """Returns a scale function that does the scaling only if
    elt_cond_fn is satisfied
    :param elt_cond_fn: function from the tuple (a, b, c, d, j, t) to a
    boolean value
    :param rest_cols_to_scale: list of indices of columns of float values
    to be scaled. These can be relative to the first float column [0, 6) or
    relative to the first column [6, 12)
    :param scalefactor: float value by which to scale rest_cols
    :return: function that takes elt and rest as arguments and returns a
    tuple of scaled "rest" values, where
    elt is the tuple (a, b, c, d, j, t) and
    rest is the tuple of associated float values,
    """
    def scaleif(elt, rest):
        if elt_cond_fn(elt):
            for col in rest_cols_to_scale:
                col_idx = col % 6
                rest[col_idx] = scalefactor * rest[col_idx]
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


def scale_off_diag_outside_valence(src, dst, nshell, scalefactor):
    scalefn = get_scaleif_fn(
        elt_cond_fn=ecf_off_diag_outside_valence(nshell),
        rest_cols_to_scale=range(6), scalefactor=scalefactor
    )
    return run(fpath_src=src, fpath_dst=dst, scalefn0=scalefn)
