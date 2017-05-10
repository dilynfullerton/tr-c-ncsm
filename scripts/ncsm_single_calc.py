#!/usr/bin/python
"""ncsm_single_calc.py

This file exists just to provide script functionality of functions in
ncsm_vce_calc.py. It is recommended that the user write their own python
script to import and use the functions in ncsm_vce_calc.py, as this would
provide greater control and clarity. However, one can use most of the
available functionality in the form of a UNIX-style command.

To run as a script (assuming this file has been added to PATH):

    $ ncsm_single_calc.py [-f] [-e] [-v] [-s SCALEFACTOR] [-t WALLTIME]
    Z A [Aeff [nmax [nshell [rm_prot [n1 n2]]]]]

In the current directory run NCSD for a single element with a single A
value and Aeff value, or (if -e) run NCSD for a range of A values with
Aeff = A.

( -f | --force )
    force recalculation of NCSD, even if
      outfiles already exist. Otherwise, calculations are not repeated.
( -e | --exact )
    A and Aeff are interpreted as Amin and Amax,
      which define a range [Amin, Amax] on which exact calculations (A=Aeff)
      are done.
( -v | --verbose )
    regular verbose output of NCSD is printed to
      stdout (unless job is submitted to cluster)
( -s | --scale-int ) SCALEFACTOR
    off-diagonal valence coupling terms of the TBME interaction are scaled
    by the value of SCALEFACTOR
( -t | --walltime ) TIME
    NCSD jobs are allotted the given walltime, a string in format hh:mm:ss
"""
from __future__ import division

from sys import argv
from collections import deque

from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException
from scripts.ncsm_vce_calc import ncsd_single_calculation
from scripts.ncsm_vce_calc import ncsd_exact_calculations
from scripts.ncsm_vce_calc import NCSD_CLUSTER_WALLTIME, NMAX, N_SHELL, N1, N2

# TODO: this script section is some really smelly code
if __name__ == "__main__":
    user_args = argv[1:]
    force0, exact_range = False, False
    progress0 = True
    cluster_submit0 = False
    walltime0 = NCSD_CLUSTER_WALLTIME
    scalefactor0 = None
    rm_prot0 = False
    nmax0 = NMAX
    nshell0 = N_SHELL
    n1_ = N1
    n2_ = N2
    while True:
        a0 = user_args[0]
        if '-f' == a0.lower() or '--force' == a0:
            force0 = True
        elif '-e' == a0 or '--exact' == a0:
            exact_range = True
        elif '-v' == a0 or '--verbose' == a0:
            progress0 = True
        elif '-s' == a0 or '--scale-int' == a0:
            user_args = user_args[1:]
            scalefactor0 = round(float(user_args[0]), 2)
        elif '-t' == a0 or '--walltime' == a0:
            user_args = user_args[1:]
            cluster_submit0 = True
            walltime0 = user_args[0]
        else:
            break
        user_args = user_args[1:]
    if not exact_range:
        if len(user_args) < 2:
            raise InvalidNumberOfArgumentsException(
                '%d' % len(user_args) +
                ' is not a valid number of arguments for ncsm_single_calc.py.' +
                ' Please enter 2 or more arguments.'
            )
        elif len(user_args) == 2:
            z0, a0 = user_args
            aeff0 = a0
        else:
            z0, a0, aeff0 = user_args[:3]
            user_args = deque(user_args[:3])
            if len(user_args) > 0:
                nmax0 = user_args.popleft()
            if len(user_args) > 0:
                nshell0 = user_args.popleft()
            if len(user_args) > 0:
                rm_prot0 = user_args.popleft()
            if len(user_args) > 1:
                n1_, n2_ = user_args[:2]
        ncsd_single_calculation(
            z=z0, a=a0, aeff=aeff0, nmax=nmax0, n1=n1_, n2=n2_, nshell=nshell0,
            scalefactor=scalefactor0, force=force0,
            cluster_submit=cluster_submit0, walltime=walltime0,
            remove_protons=rm_prot0,
        )
    else:
        if len(user_args) < 2:
            raise InvalidNumberOfArgumentsException(
                '%d' % len(user_args) +
                ' is not a valid number of arguments for ncsm_single_calc.py.' +
                ' Please enter 2 or more arguments.'
            )
        elif len(user_args) == 2:
            z0, amin = user_args
            arange0 = [amin]
        else:
            z0, amin, amax = user_args[:3]
            arange0 = range(amin, amax + 1)
            user_args = deque(user_args[:3])
            if len(user_args) > 0:
                nmax0 = user_args.popleft()
            if len(user_args) > 0:
                nshell0 = user_args.popleft()
            if len(user_args) > 0:
                rm_prot0 = user_args.popleft()
            if len(user_args) > 1:
                n1_, n2_ = user_args[:2]
        ncsd_exact_calculations(
            z=z0, a_range=arange0, nmax=nmax0, nshell=nshell0, n1=n1_, n2=n2_,
            int_scalefactor=scalefactor0, force=force0, progress=progress0,
            cluster_submit=cluster_submit0, walltime=walltime0,
            remove_protons=rm_prot0,
        )
