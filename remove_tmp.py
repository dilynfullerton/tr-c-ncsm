#!/usr/bin/python

import re
from os import walk, path, getcwd, remove


for root, dirnames, filenames in walk(getcwd()):
    for fname in filenames:
        if re.match('.*\.tmp', fname):
            fpath = path.join(root, fname)
            # print 'Removing file %s' % fpath
            remove(fpath)
            
