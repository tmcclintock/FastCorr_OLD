"""
An example of how to interface with fast_corr.

It is compared to how one might use the hankel.py
package. fastcorr is always faster
than hankel.py.

fastcorr does not benefit from parallelization
because thread creation overhead dominates the runtime.

fastcorr makes assumptions about the power
spectrum at very small and very large scales, which
allows for xi(R) to be calculated at the edges, especially
in the 1halo region. If hankel is asked to calculate xi(R)
at very small radii it will give very bad answers 
due to the assumed behavior of the scipy splines.
"""
import sys,time
sys.path.insert(0,"src/")
import fastcorr
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from hankel import SphericalHankelTransform as sht
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text',usetex=True, fontsize=20)

klin = np.genfromtxt("test_data/klin.txt")
plin = np.genfromtxt("test_data/plin.txt")
knl = np.genfromtxt("test_data/knl.txt")
pnl = np.genfromtxt("test_data/pnl.txt")

#Define the domain
NR = 1000
R = np.linspace(.1,120,NR)

N,h = 300,0.005

#Do the linear comparison
start = time.time()
xi_lin = fastcorr.calc_corr(R,klin,plin,N,h)
end = time.time()
print "fastcorr LIN time:",end-start

xi_lin2 = np.zeros_like(xi_lin)
err = np.zeros_like(xi_lin)
start = time.time()
pspl = IUS(klin,plin)
ht = sht(0,N,h)
for i in range(len(R)):
    r = R[i]
    f = lambda x: x**2*pspl(x/r)
    xi_lin2[i],err[i] = ht.transform(f)
    xi_lin2[i]/=(np.pi**2*r**3)
    err[i]/=(np.pi**2*r**3)
    continue
end = time.time()
print "hankel LIN time:",end-start

#Now do the nonlinear part
start = time.time()
xi_nl = fastcorr.calc_corr(R,knl,pnl,N,h)
end = time.time()
print "fastcorr NL time:",end-start

xi_nl2 = np.zeros_like(xi_nl)
err = np.zeros_like(xi_nl)
start = time.time()
pspl = IUS(knl,pnl)
ht = sht(0,N,h)
for i in range(len(R)):
    r = R[i]
    f = lambda x: x**2*pspl(x/r)
    xi_nl2[i],err[i] = ht.transform(f)
    xi_nl2[i]/=(np.pi**2*r**3)
    err[i]/=(np.pi**2*r**3)
    continue
end = time.time()
print "hankel NL time:",end-start

plt.loglog(R,R**2*xi_lin,label=r"$\xi_{\rm lin}$ fastcorr")
plt.loglog(R,R**2*xi_lin2,label=r"$\xi_{\rm lin}$ hankel")
plt.loglog(R,R**2*xi_nl,ls="--",label=r"$\xi_{\rm nl}$ fastcorr")
plt.loglog(R,R**2*xi_nl2,ls="--",label=r"$\xi_{\rm nl}$ hankel")
plt.xlabel(r"$R\ [{\rm Mpc}/h]$")
plt.ylabel(r"$R^2\xi(R)\ [{\rm Mpc^2}/h^2]$")
plt.ylim(.1,200)
plt.legend(loc="lower left")
plt.subplots_adjust(bottom=0.15)
plt.show()
