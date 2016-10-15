"""
An example of how to interface with fast_corr.

It is compared to how one might use the hankel.py
package. My implementation is consistently faster
for the linear power spectrum and sometimes
slightly slower for the nonlinear
power spectrum.

On top of this, it makes assumptions about the power
spectrum at very small and very large scales, which
allows for xi(R) to be sensibly calculated, especially
in the 1halo region. On the other hand,
if hankel is asked to calculate xi(R)
at very small radii it will give very bad answers,
but this is actually due to the assumed behavior
of the scipy splines.
"""
import sys,time
sys.path.insert(0,"src/")
import fast_corr
import numpy as np
import matplotlib.pyplot as plt

k = np.genfromtxt("test_data/knl.txt")
p = np.genfromtxt("test_data/pnl.txt")
NR = 500
R = np.linspace(1,120,NR)

N,h = 300,0.005

start = time.time()
xi = fast_corr.calc_corr(R,k,p,N,h)
end = time.time()
print "My time:",end-start

from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from hankel import SphericalHankelTransform as sht

xi2 = np.zeros_like(xi)
err = np.zeros_like(xi)

start = time.time()
pspl = IUS(k,p)
h = sht(0,N,h)
for i in range(len(R)):
    r = R[i]
    f = lambda x: x**2*pspl(x/r)
    fres = h._f(f,h.x)
    summation = np.pi * h.w * fres * h.j * h.dpsi
    ret = [np.sum(summation)]
    xi2[i],err[i] = h.transform(f)
    xi2[i]/=(np.pi**2*r**3)
    err[i]/=(np.pi**2*r**3)
    continue
end = time.time()
print "hankel time:",end-start

plt.loglog(R,xi)
plt.loglog(R,xi2)
plt.show()
