#!/usr/bin/python
import re
from os import walk, path, remove, getcwd

DNAME_REGEX_LIST = [
    '.*_Nhw[2-9]_',
    '.*_Nhw10_',
    'he9.*_Nhw11_',
    'he10.*_Nhw12_',
    # 'he10_.*_Nhw10_',
]
TMP_REGEX = '.*\.tmp'
RESULTS_DIR = path.join(getcwd(), 'results/ncsd')

count = 0
for root, dirs, files in walk(RESULTS_DIR):
    dname_root = path.split(root)[1]
    for rgx in DNAME_REGEX_LIST:
        if re.match(rgx, dname_root):
            break
    else:
        continue
    for fname in files:
        fpath = path.join(root, fname)
        if re.match(TMP_REGEX, fpath):
            remove(fpath)
            count += 1
print 'Removed %d files' % count
