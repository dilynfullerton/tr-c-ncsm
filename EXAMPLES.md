# tr-c-ncsm
### Examples of usage
#### Example 1: Basic _p_-shell calculation on the cluster
##### Python:

   ```python
   from ncsm_vce_calc import ncsd_vce_calculations
   from ncsm_vce_calc import ncsd_exact_calculations
   
   A_PRESCRIPTIONS = [(x,) * 3 for x in range(4, 11)] + [(4, 5, 6)]
   
   ncsd_vce_calculations(
       a_prescriptions=A_PRESCRIPTIONS, a_range=range(4, 11), nmax=6,
       cluster_submit=True, walltime='01:00:00'
   )
   ncsd_exact_calculations(
       z=2, a_range=range(7, 11), nmax=6,
       cluster_submit=True, walltime='03:00:00'
   )
   ```

The call to `ncsd_vce_calculations()` above will submit the necessary
`NCSD` jobs to the cluster using `qsub`. Then a valence cluster
expansion will be attempted on the reulsts.
Unfortunately it will fail the first
time, since the NCSD jobs will not yet be complete.
This is a limitation in submitting to the cluster, as
a method for waiting for completion has not been developed. Thus the
`ncsd_vce_calculations()` method will have to be run again once the
`NCSD` jobs have finished in order to do the `TRDENS` and VCE parts of
the calculation.

The call to `ncsd_exact_calculations()` above will submit the
remaining exact `NCSD` calculations A = Aeff = 7, 8, 9, 10 in order to
make comparisons for the full mass range.
   

##### Shell:

   ```bash
   python ncsm_vce_calc.py -t 01:00:00 -e 4 10 4 10 6;
   python ncsm_vce_calc.py -t 01:00:00 4 5 6 4 10 6;
   python ncsm_single_calc.py -t 03:00:00 -e 2 7 10 6;
   ```

_Coming soon_: Explanation

#### Example 2: Basic _sd_-shell calculation on the cluster
##### Python:

   ```python
   from ncsm_vce_calc import ncsd_vce_calculations
   from ncsm_vce_calc import ncsd_exact_calculations
   
   A_PRESCRIPTIONS = [(x,) * 3 for x in range(16, 25)] + [(16, 17, 18)]
    
   ncsd_vce_calculations(
       a_prescriptions=A_PRESCRIPTIONS, a_range=range(16, 25),
       nmax=2, nshell=2,
       cluster_submit=True, walltime='03:00:00'
   )
   ncsd_exact_calculations(
       z=8, a_range=range(19, 25), nmax=2, nshell=2,
       cluster_submit=True, walltime='05:00:00'
   )
   ```

_Coming soon_: Explanation

##### Shell:

   ```bash
   python ncsm_vce_calc.py -t 03:00:00 -e 16 24 16 24 2 2;
   python ncsm_vce_calc.py -t 03:00:00 16 17 18 16 24 2 2;
   python ncsm_single_calc.py -t 05:00:00 -e 8 19 24 2 2;
   ```

_Coming soon_: Explanation

_Coming soon_: More examples
