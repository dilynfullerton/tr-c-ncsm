#!/bin/bash
rsync -Hrl --exclue="*.tmp" --ignore-existing cougar:~/NCSM/calc/mcalc/results ./
