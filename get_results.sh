#!/bin/bash
from='cougar:~/NCSM/calc/mcalc/results';
to='./';
rsync -Hrl \
  --exclude="*.tmp" \
  --exclude="*.egv" \
  $from $to;
