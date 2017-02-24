#!/bin/bash
# puts python scripts in the scripts directory onto cougar
from='templates/*.sh';
to='cougar:~/NCSM/calc/mcalc/templates/';
rsync -r ${from} ${to};
