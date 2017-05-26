#!/usr/bin/python

from ncsm_vce_calc import ncsd_vce_calculations, ncsd_exact_calculations
from ncsm_vce_calc import generate_exact

# Script
# AP_RANGE = range(6, 11+1)
AP_RANGE = range(6, 11+1)
EX_RANGE = list(AP_RANGE)[3:]
# A_PRESCRIPTIONS = list(generate_exact(AP_RANGE, 3)) + [[6, 7, 8]]
# A_PRESCRIPTIONS = list(generate_exact(AP_RANGE, 3)) + [[4, 5, 6]]
A_PRESCRIPTIONS = list(generate_exact(AP_RANGE, 3)) + [[4, 5, 6], [6, 7, 8]]
Z = 3  # Lithium
NMAX = 2
N_SHELL = 1  # p shell
BETA_CM = 10.0
SCALEFACTOR = None
NUM_ITER = 500
WALLTIME_PRESC = '4:00:00'
WALLTIME_EXACT = '8:00:00'
USE_MPI = True

ncsd_vce_calculations(
    a_prescriptions=A_PRESCRIPTIONS,
    a_range=AP_RANGE,
    z=Z,
    nmax=NMAX,
    nshell=N_SHELL,
    ncomponent=2,
    n1=15, n2=15,
    int_scalefactor=SCALEFACTOR,
    beta_cm=BETA_CM,
    num_iter=NUM_ITER,
    cluster_submit=True,
    walltime=WALLTIME_PRESC,
    use_mpi=USE_MPI,
    remove_protons=False,
    force_ncsd=False,
    force_trdens=False,
    force_all=False,
    verbose=True,
)

ncsd_exact_calculations(
    a_range=EX_RANGE,
    z=Z,
    nmax=NMAX,
    nshell=N_SHELL,
    n1=15, n2=15,
    int_scalefactor=SCALEFACTOR,
    beta_cm=BETA_CM,
    num_iter=NUM_ITER,
    cluster_submit=True,
    walltime=WALLTIME_EXACT,
    use_mpi=USE_MPI,
    force=False,
    verbose=True,
    remove_protons=False,
)
