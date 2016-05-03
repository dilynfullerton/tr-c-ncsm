# tr-c-ncsm
Coming soon: MORE information about this stuff. Oooh weee! I can't wait

### Overview

### Dependencies
* python 2 (2.4 or higher)

### Setup
The following should be done to prepare for usage:

1. Clone this repository to the desired calculation directory.  
```bash
git clone https://github.com/dilynfullerton/tr-c-ncsm.git
```

2. Add the `NCSD` and `TRDENS` files to `PATH` in your `.bash_profile`
or `.profile`, replacing the ellipses as necessary.  
```bash
PATH=$PATH:~/.../NCSM/src/
```

3. If you want to run ncsm\_vce\_calc.py and ncsm\_single\_calc.py as
scripts, add the following to your `.bash_profile` or `.profile`, replaces
the ellipses as necessary.  
```bash
PATH=$PATH:/.../scripts/
```

4. If you want to receive emails upon completion of jobs submitted to
the cluster, add the following to `~/.../templates/job.sh`.  
```bash
#PBS -M youremail@here
```

### Usage

### Future developments
