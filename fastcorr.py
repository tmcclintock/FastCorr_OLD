"""
This is fastcorr. It calculates the matter correlation function quickly.
It can also calculate terms in RSD models.

It interfaces with the compiled fastcorr.c using ctypes.

Compile by running `python setup.py' in this directory.
"""
import numpy as np
import os, inspect
from ctypes import c_double,c_int,POINTER,cdll
#Find the library no matter where we are.
sopath = os.path.join(os.path.dirname(__file__),"_fastcorr.so")

def calc_corr(R,k,P,N=2**10,h=2**-10):
    """Calculate the 3D matter correlation function.

    Args:
        R (array_like): Radial distances; Mpc/h or Mpc.
        k (array_like): Wavenumbers; h/Mpc or Mpc^-1.
        P (array_like): Matter power spectrum; (h/Mpc)^3 or Mpc^-3.
        N (int): Number of quadrature roots; default is 2^10.
        h (float): Step size of quadrature rule; default is 2^-10.

    Returns:
        xi (array_like): Matter correlation function.

    """

    R = R.copy()
    k = k.copy()
    P = P.copy()
    cclib = cdll.LoadLibrary(sopath)
    ccc = cclib.calc_corr
    ccc.restype = c_int

    """
    Argument order for the c code:
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

def calc_xi2(R,k,P,N=2**10,h=2**-10):
    """Calculate the integral of P and j_2, for use in RSD calculations.

    Args:
        R (array_like): Radial distances; Mpc/h or Mpc.
        k (array_like): Wavenumbers; h/Mpc or Mpc^-1.
        P (array_like): Matter power spectrum; (h/Mpc)^3 or Mpc^-3.
        N (int): Number of quadrature roots; default is 2^10.
        h (float): Step size of quadrature rule; default is 2^-10.

    Returns:
        xi_2 (array_like): First perturbation for RSD.

    """
    cclib = cdll.LoadLibrary(sopath)
    ccc = cclib.calc_xi2
    ccc.restype = c_int

    """
    Argument order for c code:
    k,P,Nk,
    R,xi2,NR,
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

    xi2 = np.zeros(NR)
    xi2_in = xi2.ctypes.data_as(POINTER(c_double))

    result = ccc(k_in,P_in,Nk,R_in,xi2_in,NR,N,h)

    if result != 0:
        raise Exception("Error message recieved in fastcorr.py")
    return xi2


def calc_xi4(R,k,P,N=2**10,h=2**-10):
    """Calculate the integral of P and j_4, for use in RSD calculations.

    Args:
        R (array_like): Radial distances; Mpc/h or Mpc.
        k (array_like): Wavenumbers; h/Mpc or Mpc^-1.
        P (array_like): Matter power spectrum; (h/Mpc)^3 or Mpc^-3.
        N (int): Number of quadrature roots; default is 2^10.
        h (float): Step size of quadrature rule; default is 2^-10.

    Returns:
        xi_4 (array_like): Second perturbation for RSD.

    """
    cclib = cdll.LoadLibrary(sopath)
    ccc = cclib.calc_xi4
    ccc.restype = c_int

    """
    Argument order for c code:
    k,P,Nk,
    R,xi4,NR,
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

    xi4 = np.zeros(NR)
    xi4_in = xi4.ctypes.data_as(POINTER(c_double))

    result = ccc(k_in,P_in,Nk,R_in,xi4_in,NR,N,h)

    if result != 0:
        raise Exception("Error message recieved in fastcorr.py")
    return xi4
