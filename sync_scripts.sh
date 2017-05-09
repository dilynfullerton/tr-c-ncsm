#!/bin/bash
# puts python scripts in the scripts directory onto cougar
from='scripts/*.py';
to='cougar:/itch/dfullerton/NCSM/calc/mcalc/scripts/';
rsync -r ${from} ${to};
