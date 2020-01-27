import iminuit
import probfit
cimport cython


from libc.math cimport exp
@cython.binding(True) # IMPORTANT: this tells Cython to dump the function signature
def sigmoid(double x, double p0, double p1, double p2, double p3):
    return p0 / (1.0 + exp(-1.0 * p2 * (x - p1) ) ) + p3

@cython.binding(True)
def exponential(double x, double p0, double p1):
    return p0 * exp(x / p1)
