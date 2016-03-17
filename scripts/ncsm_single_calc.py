#!/usr/bin/python
"""ncsm_single_calc.py

To run as a script:

    $ ncsm_single_calc.py [-[Ff]] Z A [Aeff [nhw [n1 n2] | n1 n2]]

In the current directory, run NCSD for a single element with a single A
value and Aeff value.

If arguments are preceded by -f or -F, force recalculation of NCSD, even if
    outfiles already exist. Otherwise, calculations are not repeated.
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

if __name__ == "__main__":
    if '-f' == argv[1].lower():
        force0 = True
        user_args = argv[2:]
    else:
        force0 = False
        user_args = argv[1:]
    if len(user_args) == 2:
        z0, a0 = [int(x) for x in user_args]
        ncsd_single_calculation(z=z0, a=a0, aeff=a0, force=force0)
    elif len(user_args) == 3:
        z0, a0, aeff0 = [int(x) for x in user_args]
        ncsd_single_calculation(z=z0, a=a0, aeff=aeff0, force=force0)
    elif len(user_args) == 4:
        z0, a0, aeff0, nhw0 = [int(x) for x in user_args]
        ncsd_single_calculation(z=z0, a=a0, aeff=aeff0,
                                nhw=nhw0, force=force0)
    elif len(user_args) == 5:
        z0, a0, aeff0, n1_0, n2_0 = [int(x) for x in user_args]
        ncsd_single_calculation(z=z0, a=a0, aeff=aeff0,
                                n1=n1_0, n2=n2_0, force=force0)
    elif len(user_args) == 6:
        z0, a0, aeff0, nhw0, n1_0, n2_0 = [int(x) for x in user_args]
        ncsd_single_calculation(z=z0, a=a0, aeff=aeff0,
                                nhw=nhw0, n1=n1_0, n2=n2_0, force=force0)
    else:
        raise InvalidNumberOfArgumentsException(
            '%d' % len(user_args) +
            ' is not a valid number of arguments for ncsm_single_calc.py.' +
            ' Please enter 2-6 arguments.'
        )
