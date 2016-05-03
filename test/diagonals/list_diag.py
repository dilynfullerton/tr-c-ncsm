import re
from os import listdir, getcwd, path

d = getcwd()

for f in listdir(path.join(d, 'templates')):
    if re.match('TBME', f):
        f0 = path.join(d, 'templates', f)
        break

with open(f0, 'r') as fr:
    diag_lines = list()
    first_line = fr.readline()
    while True:
        line = fr.readline()
        if not line:
            break
        else:
            ldat = [int(x) for x in line.split()[:4]]
            try:
                if ldat[0] == ldat[1] == ldat[2] == ldat[3]:
                    diag_lines.append((line, ldat))
            except IndexError:
                print(line)
                raise

with open(path.join(d, 'diag.txt'), 'w') as fw:
    for line in sorted(diag_lines, key=lambda l: l[1]):
        fw.write(line[0])
