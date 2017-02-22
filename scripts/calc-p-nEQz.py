#!/usr/bin/python

from ncsm_vce_calc import ncsd_vce_calculations, ncsd_exact_calculations

AP_RANGE = range(4, 15+1)
NMAX = 2
N_SHELL = 1  # p shell

for z0 in AP_RANGE:
    ncsd_vce_calculations(
        a_prescriptions=[[z0, z0, z0]],
        a_range=[z0],
        z=z0,
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
        a_range=[z0],
        z=z0,
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
