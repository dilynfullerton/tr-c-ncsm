#!/usr/bin/python
"""ncsm_single_calc.py

To run as a script:

    $ ncsm_single_calc.py [-f] [-e] [-v] [-s] [-t walltime]
    Z A [Aeff [nhw [n1 n2 [nshell [ncomponent]]] | n1 n2]]

In the current directory run NCSD for a single element with a single A
value and Aeff value, or (if -e) run NCSD for a range of A values with
Aeff = A.

-f or -F
    force recalculation of NCSD, even if
      outfiles already exist. Otherwise, calculations are not repeated.
-e
    A and Aeff are interpreted as Amin and Amax,
      which define a range [Amin, Amax] on which exact calculations (A=Aeff)
      are done.
    nhw is interpreted as nmax
-v
    regular verbose output of NCSD is printed to
      stdout (unless job is submitted to cluster)
-s
    NCSD jobs are submitted to cluster
-t walltime
    NCSD jobs are allotted the given walltime, a string in format hh:mm:ss

If 2 arguments are given, these are Z A. (Aeff is taken to be equal to A).
If 3 arguments are given, these are Z A Aeff.
If 4 arguments are given, these are Z A Aeff nhw.
If 5 arguments are given, these are Z A Aeff n1 n2.
If 6 arguments are given, these are Z A Aeff nhw n1 n2.

If -e and...
If 7 arguments are given, these are Z A Aeff nhw n1 n2 nshell.
If 8 arguments are given, these are Z A Aeff nhw n1 n2 nshell ncomponent.
"""
from __future__ import division

from sys import argv

from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException
from scripts.ncsm_vce_calc import ncsd_single_calculation
from scripts.ncsm_vce_calc import ncsd_exact_calculations
from scripts.ncsm_vce_calc import NCSD_CLUSTER_WALLTIME

# todo this script section is some really smelly code
if __name__ == "__main__":
    user_args = argv[1:]
    force0, exact_range = False, False
    verbose0, progress0 = False, True
    cluster_submit0 = False
    walltime0 = NCSD_CLUSTER_WALLTIME
    while True:
        a0 = user_args[0]
        if '-f' == a0.lower():
            force0 = True
        elif '-e' == a0:
            exact_range = True
        elif '-v' == a0:
            verbose0, progress0 = True, False
        elif '-s' == a0:
            cluster_submit0 = True
        elif '-t' == a0:
            user_args = user_args[1:]
            cluster_submit0 = True
            walltime0 = user_args[0]
        else:
            break
        user_args = user_args[1:]
    if not exact_range:
        if len(user_args) == 2:
            z0, a0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=a0,
                force=force0, verbose=verbose0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 3:
            z0, a0, aeff0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=aeff0,
                force=force0, verbose=verbose0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 4:
            z0, a0, aeff0, nhw0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=aeff0, nhw=nhw0,
                force=force0, verbose=verbose0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 5:
            z0, a0, aeff0, n1_0, n2_0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=aeff0, n1=n1_0, n2=n2_0,
                force=force0, verbose=verbose0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 6:
            z0, a0, aeff0, nhw0, n1_0, n2_0 = [int(x) for x in user_args]
            ncsd_single_calculation(
                z=z0, a=a0, aeff=aeff0, nhw=nhw0, n1=n1_0, n2=n2_0,
                force=force0, verbose=verbose0,
                cluster_submit=cluster_submit0, walltime=walltime0,
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
                force=force0, verbose=verbose0, progress=progress0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 3:
            z0, amin, amax = [int(x) for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1),
                force=force0, verbose=verbose0, progress=progress0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 4:
            z0, amin, amax, nmax0 = [int(x) for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1), nmax=nmax0,
                force=force0, verbose=verbose0, progress=progress0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 5:
            z0, amin, amax, n1_, n2_ = [int(x) for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1), n1=n1_, n2=n2_,
                force=force0, verbose=verbose0, progress=progress0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 6:
            z0, amin, amax, nmax0, n1_, n2_ = [int(x) for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1), nmax=nmax0, n1=n1_, n2=n2_,
                force=force0, verbose=verbose0, progress=progress0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 7:
            z0, amin, amax, nmax0, n1_, n2_, nshell0 = [int(x)
                                                        for x in user_args]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1), nmax=nmax0, n1=n1_, n2=n2_,
                nshell=nshell0,
                force=force0, verbose=verbose0, progress=progress0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        elif len(user_args) == 8:
            z0, amin, amax, nmax0, n1_, n2_, nshell0, ncomp0 = [
                int(x) for x in user_args
                ]
            ncsd_exact_calculations(
                z=z0, a_range=range(amin, amax+1), nmax=nmax0, n1=n1_, n2=n2_,
                nshell=nshell0, ncomponent=ncomp0,
                force=force0, verbose=verbose0, progress=progress0,
                cluster_submit=cluster_submit0, walltime=walltime0,
            )
        else:
            raise InvalidNumberOfArgumentsException(
                '%d' % len(user_args) +
                ' is not a valid number of arguments for ncsm_single_calc.py.' +
                ' Please enter 2-6 arguments.'
            )
