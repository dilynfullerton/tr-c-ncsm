#!/usr/bin/python

from ncsm_vce_calc import ncsd_vce_calculations, ncsd_exact_calculations
from ncsm_vce_calc import _exact

# Script
AP_RANGE = range(4, 9+1)
A_PRESCRIPTIONS = _exact(AP_RANGE, 3) + [[4, 5, 6]]
Z = 2  # Helium
NMAX = 2
N_SHELL = 1  # p shell


ncsd_vce_calculations(
    a_prescriptions=A_PRESCRIPTIONS,
    a_range=AP_RANGE,
    z=Z,
    nmax=NMAX,
    nshell=N_SHELL,
    ncomponent=2,
    n1=15, n2=15,
    int_scalefactor=1.0,
    cluster_submit=True,
    walltime='48:00:00',
    remove_protons=False,
    force_ncsd=False,
    force_trdens=True,
    force_all=False,
    verbose=False,
    progress=False,
)

ncsd_exact_calculations(
    a_range=AP_RANGE,
    z=Z,
    nmax=NMAX,
    nshell=N_SHELL,
    n1=15, n2=15,
    int_scalefactor=1.0,
    cluster_submit=True,
    walltime='48:00:00',
    force=False,
    verbose=False,
    progress=False,
    remove_protons=False,
)
