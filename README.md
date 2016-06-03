# tr-c-ncsm
_Coming soon_: MORE information about this stuff. Oooh weee! I can't wait

### Overview
---
_Coming soon._

### Dependencies
---
* python 2 (2.4 or higher)

### Setup
---
The following should be done to prepare for usage:

1. Clone this repository to the desired calculation directory.  

   ```bash
   git clone https://github.com/dilynfullerton/tr-c-ncsm.git
   ```
2. Replace the `it-code-111815.f` file that is already present in the
`~/.../NCSM/src` directory with the one included here in
`./nuclear_codes` with my modifications. Recompile.  

3. Add the directory containing the `NCSD` and `TRDENS` files to
`PATH` in your `.bash_profile` or `.profile`, replacing the ellipses
as necessary.  

   ```bash
   PATH=$PATH:~/.../NCSM/src/
   ```
   
4. If you want to run `ncsm_vce_calc.py` and `ncsm_single_calc.py` as
scripts, add the following to your `.bash_profile` or `.profile`, replacing
the ellipses as necessary.  

   ```bash
   PATH=$PATH:~/.../scripts/
   ```
5. If you want to receive emails upon completion of jobs submitted to
the cluster, add the following to `./templates/job.sh`.  

   ```bash
   #PBS -m ae
   #PBS -M youremail@here
   ```
The arguments `a` and `e` prompt sending a message when the job either
(`a`) aborts or (`e`) exits normally. 

### Usage
---
#### Importing python functions
It is recommended that the user create their own python script to
import and use the primary user functions:

* `ncsd_vce_calculations()`
* `ncsd_single_calculation()`
* `ncsd_exact_calculations()`
* `ncsd_multiple_calculations()`
* `vce_single_calculation()`
* `vce_multiple_calculations()`

#### Running as a script
The user may opt to run these via a UNIX-style command. Two commands
are provided:

* `ncsm_vce_calc.py`
* `ncsm_single_calc.py`

##### Function descriptions
###### `ncsm_vce_calc.py`
   ```bash
   ncsm_vce_calc.py [-f [F_ARGS]] [-v] [-s SCALEFACTOR] [-t WALLTIME]
   ( AEFF4 AEFF5 AEFF6 | (-e|-m|-M) AP_MIN AP_MAX )
   A_MIN [A_MAX [NMAX [NSHELL [NCOMPONENT [Z [N1 N2 [RM_PROT]]]]]]]
   ```
* `(-v | --verbose)`:
    Print regular `NCSD` and `TRDENS` output to standard out.
* `(-f | --force) [FARGS]`:
    Forces recalculation of `NCSD`, `TRDENS`, and VCE even if output
    files exist. If `t` in `FARGS`, only recalculates `TRDENS` and VCE.
* `(-s | --scale-int) SCALEFACTOR`:
    Scale off-diagonal terms that couple core to valence, core to
    outside, and valence to outside terms in the TBME interaction by
	the given `SCALEFACTOR`.
* `(-t | --walltime) WALLTIME`:
    Submit the job to the cluster using `qsub`, allowing the given
    `WALLTIME`.
* `(-e | --exact) AP_MIN AP_MAX`:
    Do `NCSD`, `TRDENS`, and VCE for all prescriptions of the form
    `(Aeff, Aeff, Aeff)` where `Aeff` ranges from `AP_MIN` to `AP_MAX`
    inclusive.
* `(-m | --combinations) AP_MIN AP_MAX`:
    Do `NCSD`, `TRDENS`, and VCE for all prescriptions of the form
    `(Aeff4, Aeff5, Aeff6)`, where
    `AP_MIN <= Aeff4 < Aeff5 < Aeff6 <= AP_MAX`.
* `(-M | --multicombinations) AP_MIN AP_MAX`:
    Do `NCSD`, `TRDENS`, and VCE for all prescriptions of the form
    `(Aeff4, Aeff5, Aeff6)`, where
    `AP_MIN <= Aeff4 <= Aeff5 <= Aeff6 <= AP_MAX`.

_Coming soon_: Grammars for these shell scripts

#### Examples
See [EXAMPLES.md](EXAMPLES.md)

#### Other useful scripts
* `get_results.sh` pulls only relevant results from a remote `$from`
location to the current directory.
* `sync_scripts.sh` puts python scripts in the `scripts` directory to
a remove location, specified by `$to`.
* `vce_to_shell_sources.sh` puts VCE results (`./results/vce`) into
the NuShellX `sources` directory in a remote location, specified by `$to`.

### Future developments
---
_Coming soon._  
