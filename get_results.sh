#!/bin/bash
# pulls relevant results from cougar to the current directory

# cougar -> here
#from='cougar:~/NCSM/calc/mcalc/results';
from='cougar:~/NCSM/calc/mcalc/try';
to='./';

# here -> itheory
# from='./results';
# to='itheory:~/workspace/tr-c-ncsm/';

rsync -Hrlv \
  --exclude="*.tmp" \
  --exclude="*.egv" \
  $from $to;
