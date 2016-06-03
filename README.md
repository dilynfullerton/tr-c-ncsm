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
Perform `NCSD`, `TRDENS`, and VCE calculations for the desired
prescriptions.  
   ```bash
   ncsm_vce_calc.py [-f [F_ARGS]] [-v] [-s SCALEFACTOR] [-t WALLTIME]
   ( AEFF4 AEFF5 AEFF6 | (-e|-m|-M) AP_MIN AP_MAX )
   A_MIN [A_MAX [NMAX [NSHELL [NCOMPONENT [Z [N1 N2 [RM_PROT]]]]]]]
   ```

* `F_ARGS ::= string`  
Character arguments for what parts of the calculation to
redo. Currently the only meaningful argument is `t`.
  * `t`: Force redo `TRDENS` and VCE parts
* `SCALEFACTOR ::= float`  
Factor by which to scale off-diagonal coupling terms in the TBME
interaction.
* `WALLTIME ::= "[0-9]{2,}:[0-9]{2}:[0-9]{2}"`  
Amount of time to give to submitted `NCSD` jobs in the format
`hh:mm:ss`, with more `h` as necessary.
* `AEFF4 AEFF5 AEFF6 ::= int int int`  
If none of `-e`, `-m`, or `-M` are used, these three integers specify
a *single* A-prescription for which to do the calculations.
* `AP_MIN AP_MAX ::= int int`  
If `-e`, `-m`, or `-M` are given, these integers specify the inclusive
domain `[AP_MIN, AP_MAX]` over which to get A-prescriptions.
* `A_MIN ::= int`  
Specifies the minimum mass number for which NuShellX \*.int files
are generated. If `A_MAX` is not specified, this will be the only mass
number for which a NuShellX \*.int file is generated.
* `A_MAX ::= int`  
With `A_MIN`, this defines an inclusive mass range `[A_MIN, A_MAX]`
for which NuShellX \*.int files are generated.
**Note:** These will all be the *same* interaction, linked to have
different file names to facilitate the use of my scripts in
[tr-c-nushellx](https://github.com/dilynfullerton/tr-c-nushellx).
* `NMAX ::= int`  
Major oscillator shell truncation level, defined by the number of
shells *more* than the minimum required.
* `NSHELL ::= int`  
Shell model nuclear shell number. (0=s, 1=p, 2=sd, ...).
* `NCOMPONENT ::= (1|2)`  
Number of components. (1 -> neutrons only, 2 -> protons & neutrons)
* `Z ::= int`  
Proton number, if different from A0/2, where A0 is the minimum mass
number in the shell.
* `N1 N2 ::= int int`  
Truncation levels for TBME interaction file.
* `RM_PROT ::= (0|1)`  
If `0` (false), proceed as usual. If `1` (true), remove the proton
parts of the interaction by scaling them to 0.

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
