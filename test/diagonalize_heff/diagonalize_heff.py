from __future__ import division, print_function
from sys import argv
import numpy as np
from numpy.linalg import eigh, eig

F_READ = 'Heff_OLS.dat'
F_WRITE = 'Heff_Eigen.dat'

def get_heff_matrix(fpath=F_READ):
    """Reads the Heff matrix from the given filepath and returns
    a 2D numpy array containing the matrix
    """
    with open(fpath, 'r') as fr:
        line = fr.readline()
        dim = int(line.split()[0])
        for i in range(dim):
            fr.readline()
        matrix = np.empty(shape=(dim, dim))
        for i in range(dim):
            ldat = fr.readline().split()
            row = np.array([float(x) for x in ldat])
            matrix[i, :] = row
    return matrix


def write_eigens(eigenvalues, eigenvectors, fpath=F_WRITE):
    n = len(eigenvalues)
    with open(fpath, 'w') as fw:
        fw.write('Eigenvalues\n')
        for eval, i in zip(eigenvalues, range(n)):
            fw.write('  {:4}  {:16.8f}\n'.format(i+1, eval))
        fw.write('\nEigenvectors\n')
        for evect, i in zip(eigenvectors.T, range(n)):
            fw.write('  {:4}\n'.format(i+1))
            for val in evect:
                fw.write('        {:16.8f}\n'.format(val))


if __name__ == '__main__':
    if len(argv) < 2:
        mat = get_heff_matrix(F_READ)
        write_eigens(*eigh(mat))
    elif len(argv) == 2:
        f_read = argv[1]
        mat = get_heff_matrix(fpath=f_read)
        write_eigens(*eigh(mat))
    elif len(argv) == 3:
        f_read, f_write = argv[1:]
        mat = get_heff_matrix(fpath=f_read)
        evals, evects = eigh(mat)
        write_eigens(
            eigenvalues=evals, eigenvectors=evects, fpath=f_write)
    else:
        print('Invalid number of arguments\n')
        
        
