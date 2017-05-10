#!/usr/bin/python

from ncsm_vce_calc import ncsd_vce_calculations, ncsd_exact_calculations

NMAX = 0
N_SHELL = 1  # p shell

for z0 in range(2, 7):
    ncsd_vce_calculations(
        a_prescriptions=[[2*z0, 2*z0, 2*z0]],
        a_range=[2*z0],
        z=z0,
        nmax=NMAX,
        nshell=N_SHELL,
        ncomponent=2,
        n1=15, n2=15,
        int_scalefactor=None,
        cluster_submit=True,
        walltime='6:00:00',
        remove_protons=False,
        force_ncsd=False,
        force_trdens=True,
        force_all=False,
        progress=False,
    )

    # ncsd_exact_calculations(
    #     a_range=[2*z0],
    #     z=z0,
    #     nmax=NMAX,
    #     nshell=N_SHELL,
    #     n1=15, n2=15,
    #     int_scalefactor=None,
    #     cluster_submit=True,
    #     walltime='12:00:00',
    #     force=False,
    #     progress=False,
    #     remove_protons=False,
    # )
