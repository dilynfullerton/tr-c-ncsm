#!/bin/bash
# puts vce results into the NuShellX sources directory
from='./results/vce';
to='~/nushellx/linux/calculations/t0/sources/';
#to='cougar:~/nushellx/linux/calculations/t0/sources/';
rsync -r ${from} ${to};
