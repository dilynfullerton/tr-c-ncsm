#!/bin/bash
rsync -Hrl --exclude="*.tmp" cougar:~/NCSM/calc/mcalc/results ./
