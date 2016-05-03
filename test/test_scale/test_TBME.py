#!/usr/bin/python
import re
from sys import argv
from os import listdir, getcwd, path


def get_tbme_in_cwd():
    for fname in listdir(getcwd()):
        if re.match('TBME', fname):
            return path.join(getcwd(), fname)
    else:
        print('TBME file not found')
        return None


def test_tbme(fpath_tbme):
    fin = open(fpath_tbme, 'r')
    fout = open(path.join(getcwd(), 'tbme_test_result'), 'w')
    fin.readline()
    for line in fin:
        ldat = line.split()
        a, b, c, d = [int(x) for x in ldat[:4]]
        scaled_cols = [float(x) for x in [ldat[6]] + ldat[-3:]]
        for col in scaled_cols:
            if abs(col) != 0.0:
                all_zero = False
                break
        else:
            all_zero = True
        # write all potentially problematic lines
        if (b == 1 < d) or (d == 1 < b):  # core to outside
            if not all_zero:  # not scaled
                fout.write(line)
        elif (b <= 3 < d) or (d <= 3 < b):  # valence to outside
            if not all_zero:  # not scaled
                fout.write(line)
        elif all_zero:  # diagonal term that is scaled
            fout.write(line)
    fout.close()
    fin.close()

if __name__ == "__main__":
    if len(argv) > 1:
        fpath_tbme0 = argv[1]
    else:
        fpath_tbme0 = get_tbme_in_cwd()
    if fpath_tbme0 is not None:
        test_tbme(fpath_tbme0)
