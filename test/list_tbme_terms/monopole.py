from os import getcwd, path

idx = (1, 3, 1, 3)
dpath = getcwd()
fname = 'nmax0'
fpath = path.join(dpath, fname)

numerator = 0.0
denominator = 0.0

with open(fpath, 'r') as f:
    for line in f:
        ldat = line.split()
        a, b, c, d, j, t = [int(x) for x in ldat[:6]]
        if (a, b, c, d) == idx:
            trel = float(ldat[6])
            a_jt = (2 * j + 1) * (2 * t + 1)
            numerator += a_jt * trel
            denominator += a_jt

print('Monopole for {}: {}'.format(idx, numerator/denominator))
        
        
