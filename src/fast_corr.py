"""
This is a python wrapper that calls the calc_corr() c function.
This interfaces through c_types so that the user
doesn't have to.
"""
import numpy as np
import ctypes
from ctypes import c_double,c_int,POINTER,cdll

def calc_corr(R,k,P,N=200,h=0.005):
    cclib = cdll.LoadLibrary("src/c_fast_corr.so")
    ccc = cclib.calc_corr
    ccc.restype = c_int

    """
    Arguments are:
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
        raise Exception("Error message recieved in fast_corr.py")
    return xi
