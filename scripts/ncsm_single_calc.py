#!/usr/bin/python
"""ncsm_single_calc.py

To run as a script:

    $ ncsm_single_calc.py [-f] [-e] Z A [Aeff [nhw [n1 n2] | n1 n2]]

In the current directory, run NCSD for a single element with a single A
value and Aeff value.

If arguments are preceded by -f or -F, force recalculation of NCSD, even if
    outfiles already exist. Otherwise, calculations are not repeated.
If arguments are preceded by -e, A and Aeff are interpreted as Amin and Amax,
    which define a range [Amin, Amax] on which exact calculations (A=Aeff)
    are done.
If 2 arguments are given, these are Z A. Aeff is taken to be equal to A.
If 3 arguments are given, these are Z A Aeff.
If 4 arguments are given, these are Z A Aeff nhw.
If 5 arguments are given, these are Z A Aeff n1 n2.
If 6 arguments are given, these are Z A Aeff nhw n1 n2.
"""
from __future__ import division

from sys import argv

from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException
from scripts.ncsm_vce_calc import ncsd_single_calculation
from scripts.ncsm_vce_calc import ncsd_exact_calculations

# todo this script section is some really smell code
if __name__ == "__main__":
    user_args = argv[1:]
    force0, exact_range = False, False
    verbose = False
    while True:
        a0 = user_args[0]
        if '-f' == a0.lower():
            force0 = True
        elif '-e' == a0.lower():
            exact_range = True
        elif '-v' == a0.lower():
            verbose = True
        else:
            break
        user_args = user_args[1:]
    if not exact_range:
        if len(user_args) == 2:
            z0, a0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=a0,
                force=force0, verbose=verbose
            )
        elif len(user_args) == 3:
            z0, a0, aeff0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=aeff0,
                force=force0, verbose=verbose
            )
        elif len(user_args) == 4:
            z0, a0, aeff0, nhw0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=aeff0, nhw=nhw0,
                force=force0, verbose=verbose
            )
        elif len(user_args) == 5:
            z0, a0, aeff0, n1_0, n2_0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=aeff0, n1=n1_0, n2=n2_0,
                force=force0, verbose=verbose
            )
        elif len(user_args) == 6:
            z0, a0, aeff0, nhw0, n1_0, n2_0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=aeff0, nhw=nhw0, n1=n1_0, n2=n2_0,
                force=force0, verbose=verbose
            )
        else:
            raise InvalidNumberOfArgumentsException(
                '%d' % len(user_args) +
                ' is not a valid number of arguments for ncsm_single_calc.py.' +
                ' Please enter 2-6 arguments.'
            )
    else:
        if len(user_args) == 2:
            z0, amin = [int(x) for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=[amin],
                force=force0, verbose=verbose
            )
        elif len(user_args) == 3:
            z0, amin, amax = [int(x) for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1),
                force=force0, verbose=verbose
            )
        elif len(user_args) == 4:
            z0, amin, amax, nhw0 = [int(x) for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1), nhw=nhw0,
                force=force0, verbose=verbose
            )
        elif len(user_args) == 5:
            z0, amin, amax, n1_, n2_ = [int(x) for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1), n1=n1_, n2=n2_,
                force=force0, verbose=verbose
            )
        elif len(user_args) == 6:
            z0, amin, amax, nhw0, n1_, n2_ = [int(x) for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1), nhw=nhw0, n1=n1_, n2=n2_,
                force=force0, verbose=verbose
            )
        else:
            raise InvalidNumberOfArgumentsException(
                '%d' % len(user_args) +
                ' is not a valid number of arguments for ncsm_single_calc.py.' +
                ' Please enter 2-6 arguments.'
            )
