"""
This is a python wrapper that calls the calc_corr() c function.
This interfaces through c_types so that the user
doesn't have to.
"""
import numpy as np
import ctypes
from ctypes import c_double,c_int,POINTER,cdll

"""
Calculate the correlation function.

R: tangential radial distance - units in either Mpc/h or Mpc
k: wavenumber - units in either h/Mpc or Mpc^-1
P: power spectrum - units in either (h/Mpc)^3 or Mpc^-3
N: number of roots of j_0 to evaluate
h: step size for the quadrature routine
"""
def calc_corr(R,k,P,N=300,h=0.005):
    cclib = cdll.LoadLibrary("src/c_fastcorr.so")
    ccc = cclib.calc_corr
    ccc.restype = c_int

    """
    Argument order:
    k,P,Nk,
    R,xi,NR,
    N,h
    """
    Nk = len(k)
    if Nk != len(P):
        raise Exception("len(k)!=len(P)")
    NR = len(R)
    ccc.argtypes=[POINTER(c_double),POINTER(c_double),c_int,\
                  POINTER(c_double),POINTER(c_double),c_int,\
                  c_int,c_double]
    k_in = k.ctypes.data_as(POINTER(c_double))
    P_in = P.ctypes.data_as(POINTER(c_double))
    R_in = R.ctypes.data_as(POINTER(c_double))

    xi = np.zeros(NR)
    xi_in = xi.ctypes.data_as(POINTER(c_double))

    result = ccc(k_in,P_in,Nk,R_in,xi_in,NR,N,h)

    if result != 0:
        raise Exception("Error message recieved in fastcorr.py")
    #Return the correlation function. It is of length len(R)
    return xi
