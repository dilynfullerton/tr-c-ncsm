# tr-c-ncsm
### Script grammars
#### ncsm_vce_calc.py
Perform `NCSD`, `TRDENS`, and VCE calculations for the desired
prescriptions.  
##### Form:
   ```bash
   ncsm_vce_calc.py [-f [F_ARGS]] [-v] [-s SCALEFACTOR] [-t WALLTIME]
   ( AEFF4 AEFF5 AEFF6 | (-e|-m|-M) AP_MIN AP_MAX )
   A_MIN [A_MAX [NMAX [NSHELL [NCOMPONENT [Z [N1 N2 [RM_PROT]]]]]]]
   ```

##### Flags:
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

##### Defintions:
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

#### ncsm_single_calc.py
Perform `NCSD` calculations only for a single mass number or for a range.
##### Form:
   ```bash
   ncsm_single_calc.py [-f] [-e] [-v] [-s SCALEFACTOR] [-t WALLTIME]
   Z A [AEFF [NMAX [NSHELL [RM_PROT [N1 N2]]]]]
   ```

##### Flags:
* `(-e | --exact)`:  
If given, `A` and `AEFF` are interpreted as `A_MIN` and `A_MAX`, which
define an inclusive range `[A_MIN, A_MAX]` on which `NCSD`
calculations are done with `AEFF = A`.
* Other flag behaviors same as above.

##### Definitions:
* `A ::= int`  
Mass number, or `A_MIN` if `-e`.
* `AEFF ::= int`  
Effective mass number, or `A_MAX` if `-e`.
* Other constants same as above.
