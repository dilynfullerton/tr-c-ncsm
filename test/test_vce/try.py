import numpy as np
from scipy.optimize import leastsq


def lsfn(params, fitfn, xarr, yarr, spe1arr, spe3arr):
    yfit = np.empty_like(yarr)
    for x, spe1, spe3, i in zip(xarr, spe1arr, spe3arr, range(len(xarr))):
        yfit[i] = fitfn(x, spe1, spe3, params)
    return yarr - yfit


def ffn(x, spe1, spe3, params):
    a1, a2, a3 = params
    # return x*int(a1)/2 + 39.757*int(a2)/2 + (spe1+spe3)*int(a3)/2
    # return x*(a1)/2 + 39.757*(a2)/2 + (spe1+spe3)*(a3)/2
    return x + 39.757 * a2 / 2 - (spe1 + spe3) * a3 / 2


X = np.array([
    0.1288338 * 10 ** 2,
    -0.7902375 * 10,
    0.2849067 * 10,
    -0.4319996,
    -0.7124451,
])

Y = np.zeros(5)

S1 = 12.9671
S3 = 4.6039

SPE1 = np.array([S1, S3, S1, S1, S3])
SPE3 = np.array([S1] + 4 * [S3])

ans = leastsq(
    func=lsfn,
    x0=np.ones(3),
    args=(ffn, X, Y, SPE1, SPE3),
    full_output=True,
)

print(ans)
params0 = ans[0]

YFIT = np.empty_like(Y)
for x_, spe1_, spe3_, i_ in zip(X, SPE1, SPE3, range(len(X))):
    YFIT[i_] = ffn(x_, spe1_, spe3_, params0)
print(YFIT + X)
