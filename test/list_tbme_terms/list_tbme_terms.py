import re
from os import listdir, getcwd, path

dpath = getcwd()
for fname in listdir(dpath):
    if re.match('TBME', fname):
        break
else:
    print('TBME file not found')
fpath = path.join(dpath, fname)
with open(fpath, 'r') as fin:
    with open(path.join(dpath, 'nmax0'), 'w') as fout:
        first_line = fin.readline()
        write_map = dict()
        for line in fin:
            ldat = line.split()
            a, b, c, d, j, t = [int(x) for x in ldat[:6]]
            space = [1, 2, 3]
            core = [1]
            val = [2, 3]
            in_space = b in space and d in space
            if not in_space:
                continue
            # if b in core and d in val:
            #     continue
            # elif not (b in core and d in val):
            #     continue
            write_map[(j, t, a, b, c, d)] = line
        for k, v in sorted(write_map.items()):
            fout.write(v)
