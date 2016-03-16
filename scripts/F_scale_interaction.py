#!/usr/bin/python
"""F_scale_interaction.py

To run as a script:

    $ F_scale_interaction.py [-[f|F]] src dst [scalefactor]

Reads the interaction TBME file from the given src file path, and writes it to
dst, changing the float values in accordance with scalefn. As it is not
feasible to pass a function when calling a script, there is not really any
point in running this as a script, except to ensure proper function.
"""

from os import path
from sys import argv
from ncsm_vce_calc import InvalidNumberOfArgumentsException

SCALE_FACTOR = 1.0
_LINE_FMT = ' %2d' * 6 + ' %11.6f' * 6 + '\n'


def get_scaleif_fn(elt_cond_fn, rest_cols, scalefactor):
    def scaleif(elt, rest):
        if elt_cond_fn(elt):
            for col in rest_cols:
                rest[col] = scalefactor * rest[col]
        return rest
    return scaleif


def run(fpath_src, fpath_dst,
        scalefn=lambda a, b: b,
        force=False,
        _line_fmt=_LINE_FMT):
    if force or not path.exists(fpath_dst):
        fin = open(fpath_src, 'r')
        fout = open(fpath_dst, 'w')
        first_line = fin.readline()
        fout.write(first_line)
        for line in fin:
            ldat = line.split()
            elt = [int(x) for x in ldat[:6]]  # a, b, c, d, j, t
            rest = [float(x) for x in ldat[6:]]
            next_rest = scalefn(elt, rest)
            next_line = _line_fmt % (tuple(elt) + tuple(next_rest))
            fout.write(next_line)
        fout.close()
        fin.close()


if __name__ == "__main__":
    user_args = argv[1:]
    if '-f' == user_args[0].lower():
        force0 = True
        user_args = user_args[1:]
    else:
        force0 = False
    if len(user_args) == 2:
        src, dst = user_args[:2]
        scalefn = get_scaleif_fn(
            elt_cond_fn=lambda elt: elt[5] == 1,
            rest_cols=[3],
            scalefactor=SCALE_FACTOR
        )
        run(src, dst, scalefn, force=force0)
    elif len(user_args) == 3:
        src, dst = user_args[:2]
        scalefn = get_scaleif_fn(
            elt_cond_fn=lambda elt: elt[5] == 1,
            rest_cols=[3],
            scalefactor=float(user_args[2])
        )
        run(src, dst, scalefn, force=force0)
    else:
        raise InvalidNumberOfArgumentsException(
            'F_scale_interaction.py requires 2 or 3 arguments but was ' +
            'provided with %d. ' % len(user_args)
        )
