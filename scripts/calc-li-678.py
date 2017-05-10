#!/usr/bin/python

from ncsm_vce_calc import ncsd_vce_calculations, ncsd_exact_calculations
from ncsm_vce_calc import generate_exact

# Script
AP_RANGE = range(4, 11+1)
EX_RANGE = list(AP_RANGE)[3:]
A_PRESCRIPTIONS = list(generate_exact(AP_RANGE, 3)) + [[6, 7, 8]]
Z = 3  # Lithium
NMAX = 0
N_SHELL = 1  # p shell

ncsd_vce_calculations(
    a_prescriptions=A_PRESCRIPTIONS,
    a_range=AP_RANGE,
    z=Z,
    nmax=NMAX,
    nshell=N_SHELL,
    ncomponent=2,
    n1=15, n2=15,
    int_scalefactor=None,
    cluster_submit=True,
    walltime='4:00:00',
    remove_protons=False,
    force_ncsd=False,
    force_trdens=False,
    force_all=False,
    progress=True,
)

ncsd_exact_calculations(
    a_range=EX_RANGE,
    z=Z,
    nmax=NMAX,
    nshell=N_SHELL,
    n1=15, n2=15,
    int_scalefactor=None,
    cluster_submit=True,
    walltime='6:00:00',
    force=False,
    progress=True,
    remove_protons=False,
)
