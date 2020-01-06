import cython
import numpy as np
cimport cython
cimport numpy as cnp
import scipy.signal as signal

@cython.boundscheck(False)
def smooth(x, window_len=11, window="hanning"):
    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
        raise ValueError(
            "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        )

    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-2 : -window_len - 1 : -1]]
    if window == "flat":  # moving average
        w = np.ones(window_len)
    else:
        w = eval("np." + window + "(window_len)")

    y = signal.convolve(s, w / w.sum(), mode="same")

    return y



@cython.boundscheck(False)
cpdef smooth_v2(cnp.ndarray x, int window_len=11, window="hanning"):
    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
        raise ValueError(
            "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        )
    cdef cnp.ndarray s = np.r_[x[window_len - 1 : 0 : -1], x, x[-2 : -window_len - 1 : -1]]
    if window == "flat":  # moving average
        w = np.ones(window_len)
    else:
        w = eval("np." + window + "(window_len)")

    cdef cnp.ndarray y = signal.convolve(s, w / w.sum(), mode="same")

    return y