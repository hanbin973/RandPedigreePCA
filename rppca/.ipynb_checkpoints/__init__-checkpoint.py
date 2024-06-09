import numpy as np
import scipy
import scipy.sparse as sparse
import numba

@numba.njit("void(i4[:], i4[:], f8[:], f8[:,:])")
def back_2dsolve_c(Lp, Li, Lx, Y):
    """
    `L` is lower-triangular Cholesky factor in CSC format: solve `L X = Y`.
    `Y` is updated in-place.
    """
    X = Y
    for j in range(0, X.shape[0]):
        for k in range(X.shape[1]):
            X[j,k] /= Lx[Lp[j]]
            for p in range(Lp[j] + 1, Lp[j + 1]):
                X[Li[p],k] -= Lx[p] * X[j,k]

@numba.njit("void(i4[:], i4[:], f8[:], f8[:,:])")
def forward_2dsolve_c(Lp, Li, Lx, Y):
    """
    `L` is lower-triangular Cholesky factor in CSC format: solve `L' X = Y`.
    `Y` is updated in-place.
    """
    X = Y
    for j in range(X.shape[0] - 1, -1, -1):
        for k in range(X.shape[1]):
            for p in range(Lp[j] + 1, Lp[j + 1]):
                X[j,k] -= Lx[p] * X[Li[p],k]
            X[j,k] /= Lx[Lp[j]]

@numba.njit("void(i4[:], i4[:], f8[:], f8[:,:])")
def back_2dsolve_f(Lp, Li, Lx, Y):
    """
    `L` is lower-triangular Cholesky factor in CSC format: solve `L X = Y`.
    `Y` is updated in-place.
    """
    X = Y
    for k in range(X.shape[1]):
        for j in range(0, X.shape[0]):
            X[j,k] /= Lx[Lp[j]]
            for p in range(Lp[j] + 1, Lp[j + 1]):
                X[Li[p],k] -= Lx[p] * X[j,k]

@numba.njit("void(i4[:], i4[:], f8[:], f8[:,:])")
def forward_2dsolve_f(Lp, Li, Lx, Y):
    """
    `L` is lower-triangular Cholesky factor in CSC format: solve `L' X = Y`.
    `Y` is updated in-place.
    """
    X = Y
    for k in range(X.shape[1]):
        for j in range(X.shape[0] - 1, -1, -1):
            for p in range(Lp[j] + 1, Lp[j + 1]):
                X[j,k] -= Lx[p] * X[Li[p],k]
            X[j,k] /= Lx[Lp[j]]

@numba.njit("Tuple((f8[:,:], f8[:,:]))(f8[:,:])")
def qr(mat):
    return np.linalg.qr(mat)

class Design:
    def __init__(self, L):
        """
        L is a CSC lower triangular matrix
        """
        self.Lp = L.indptr
        self.Li = L.indices
        self.Lx = L.data

        self.num_inds = L.shape[0]
    
    def dot_left(self, left):
        if left.flags['C_CONTIGUOUS']:
            back_2dsolve_c(self.Lp, self.Li, self.Lx, left)
        elif left.flags['F_CONTIGUOUS']:
            back_2dsolve_f(self.Lp, self.Li, self.Lx, left)
        else:
            raise ValueError('left is not contiguous')        
        return left.copy()

    def dot_right(self, right):
        out = right.copy()
        if out.flags['C_CONTIGUOUS']:
            forward_2dsolve_c(self.Lp, self.Li, self.Lx, out)
        elif out.flags['F_CONTIGHOUS']:
            forward_2dsolve_f(self.Lp, self.Li, self.Lx, out)
        else:
            raise ValueError('out is not contiguous')
        return out

def randomized_svd(design, n_components=2, n_iter=5, n_oversamples=5, random_matrix=None, seed=None):
    """
    design is the RowEdgeDesign class object
    """
    
    if random_matrix is None:
        if seed is None:
            seed = 0
        np.random.seed(seed)
        random_matrix = np.random.normal(size=(n_components+n_oversamples, design.num_inds)).T

    sample_matrix = design.dot_left(random_matrix)
    range_old, _ = qr(sample_matrix)
    
    for i in range(n_iter):
        range_new, _ = qr(design.dot_right(range_old))
        range_old, _ = qr(design.dot_left(range_new))
        
    U, S, V = np.linalg.svd(design.dot_right(range_old).T, full_matrices=False)
    return (range_old @ U)[:,:n_components], S[:n_components], V[:n_components,:]