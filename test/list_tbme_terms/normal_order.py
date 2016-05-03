from os import getcwd, path
from math import sqrt


def get_n(a):
    return int(a == 1)


def get_j(i):
    if i == 1:
        return 0
    elif i == 2:
        return 1/2
    elif i == 3:
        return 3/2


def normal_order(v_map, A):
    e0 = 0
    f_ij = dict()
    g_ijkl = dict()
    for label, v in sorted(v_map.items()):
        a, b, c, d, j, t = label
        norm = sqrt((1+int(a==b))*(1+int(c==d))) 
        v = v * norm * 2 / (A - 1)
        if (a, b) == (c, d):
            if get_n(a) * get_n(b):
                print('{} {:8.4f}'.format(label, v))
            e0 += 1/2 * get_n(a) * get_n(b) * (2*j+1) * (2*t+1) * v
        if a == c:
            if (b, d) not in f_ij:
                f_ij[(b, d)] = 0
            f_ij[(b, d)] += get_n(a) * (2*j+1) * (2*t+1) / (2*get_j(b)+1) * v
        g_ijkl[(a, b, c, d, j, t)] = v
    print()
    return e0, f_ij, g_ijkl


def write_file(fpath, e0, f_ij, g_ijkl):
    with open(fpath, 'w') as fw:
        fw.write('Zero body\n')
        fw.write('               {:8.4f}\n'.format(e0))
        fw.write('One body\n')
        for ij, f in sorted(f_ij.items()):
            if f == 0:
                continue
            i, j = ij
            fw.write('  {} {}:         {:8.4f}\n'.format(i, j, f))
        fw.write('Two body\n')
        for idx, g in sorted(g_ijkl.items()):
            if g == 0:
                continue
            i, j, k, l, J, T = idx
            fw.write('  {} {} {} {} {} {}: {:8.4f}\n'
                     ''.format(i, j, k, l, J, T, g))


dpath = getcwd()
fname = 'nmax0'
fpath = path.join(dpath, fname)
trel_map = dict()
vpn_map = dict()
vpp_map = dict()
vnn_map = dict()

# fill maps
with open(fpath, 'r') as f:
    for line in f:
        ldat = line.split()
        label = tuple([int(x) for x in ldat[:6]])
        trel, hrel, clmb, vpn, vpp, vnn = [float(x) for x in ldat[6:]]
        trel_map[label] = trel
        vpn_map[label] = vpn
        vpp_map[label] = vpp
        vnn_map[label] = vnn

A = 4
Aeff = 4

trel_e0, trel_f_ij, trel_g_ijkl = normal_order(trel_map, A)
vnn_e0, vnn_f_ij, vnn_g_ijkl = normal_order(vnn_map, A)
# vpn_e0, vpn_f_ij, vpn_g_ijkl = normal_order(vpn_map)
# vpp_e0, vpp_f_ij, vpp_g_ijkl = normal_order(vpp_map)

write_file(path.join(dpath, 'NO_trel'), trel_e0, trel_f_ij, trel_g_ijkl)
write_file(path.join(dpath, 'NO_vnn'), vnn_e0, vnn_f_ij, vnn_g_ijkl)

print('Trel   = {:8.4f}'.format(trel_e0))
print('Vnn    = {:8.4f}'.format(vnn_e0))
print()
he4_4 = (1 - 1/A) * trel_e0 + vnn_e0
print('A = {}, Aeff = {}'.format(A, Aeff))
print('  Energy = {:8.4f}'.format(he4_4))
