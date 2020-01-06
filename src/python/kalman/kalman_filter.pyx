from cython.parallel import prange
import cython
import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False) 
cpdef kalman_filter(np.ndarray[np.float32_t, ndim=1] z, Py_ssize_t n_iter=10):
    """
    Variables
    ===
    array of the input np.ndarray[np.float32_t, ndim=1] 
    n_iter = number of iterations (default 10)
    
    """
#   DTYPE =
    cdef Py_ssize_t sz = z.shape[0]
    Q = 1e-5 # process variance
    temp = np.zeros(sz, dtype=np.float32)
    # allocate space for arrays
    cdef float[:] xhat = temp                                    # a posteri estimate of x
    cdef float[:] P = np.zeros(sz, dtype=np.float32)           # a posteri error estimate
    cdef float[:] xhatminus = np.zeros(sz, dtype=np.float32)   # a priori estimate of x
    cdef float[:] Pminus = np.zeros(sz, dtype=np.float32)      # a priori error estimate
    cdef float[:] K = np.zeros(sz, dtype=np.float32)           # gain or blending factor

    R = 0.01**2  # estimate of measurement variance, change to see effect

    # intial guesses
    xhat[0] = 0.0
    P[0] = 2.0
    cdef Py_ssize_t k, n_it = n_iter
    for k in prange(1, n_it, nogil=True):
        # time update
        xhatminus[k] = xhat[k-1]
        Pminus[k] = P[k-1]+Q

        # measurement update
        K[k] = Pminus[k] / ( Pminus[k] + R )
        xhat[k] = xhatminus[k] + K[k] * ( z[k] - xhatminus[k] )
        P[k] = ( 1 - K[k] ) * Pminus[k]
    
    return temp
