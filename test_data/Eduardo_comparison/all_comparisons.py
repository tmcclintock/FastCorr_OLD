"""
This script makes all of the relevant comparisons between
the various P(k) and Xi(R) curves supplied by
Eduardo Rozo and Matt Becker.

Each comparison consists of:
1) Load in P(k)
2) Look at plot of P(k)
3) Compute Xi(R)
4) Look at plot of Xi(R) comparison

NOTE: Eduardo has actually supplied Delta^2 = P^2 k^3/ 2pi^2
"""
import fastcorr
import numpy as np
import time
import matplotlib.pyplot as plt

N=32000
h=0.000006125
N=2**12
h=2**-10#-15
print "Product = %f"%(N*h)
print "\tN=%d\n\th=%f"%(N,h)

def call_fc(R,k,P):
    return fastcorr.calc_corr(R,k,P,N,h)

def comparison(pname,xiname,name="given",see_P=False):
    k,Del = np.genfromtxt(pname).T
    P = Del*2*np.pi**2/k**3
    if see_P:
        plt.loglog(k,P,label="Theirs")
        plt.show()
        plt.clf()
    R,R2xi = np.genfromtxt(xiname).T
    xi = R2xi/R**2
    start = time.time()
    xi_fc = call_fc(R,k,P)
    end = time.time()
    print "Time in sections = %.5f"%(end - start)
    pdiff = (xi-xi_fc)/xi * 100
    f,axes = plt.subplots(2,sharex=True)
    axes[0].loglog(R,xi,label=name)
    axes[0].loglog(R,xi_fc,'b--',label="FastCorr")
    axes[0].legend(loc='lower left')
    axes[0].set_xlim(min(R),100)
    axes[1].plot(R,pdiff)
    lim = 10
    axes[1].plot([min(R),100],[1,1],'k:')
    axes[1].plot([min(R),100],[-1,-1],'k:')
    axes[1].set_ylim(-lim,lim)
    axes[1].set_ylabel("% Diff")
    plt.show()
    plt.clf()
    plt.close()
    return
"""
1) Comparison with Eduardo's Xi_lin at z=0
2) Comparison with CAMB's Xi_lin at z=0
3) Comparison with Eduardo's Xi_lin at z=1
4) Comparison with Eduardo's Xi_nonlin at z=0
5) Comparison with Eduardo's Xi_nonlin at z=1
"""
comparison("ps_z0.0.dat","r2xilin_z0.0.dat","Eduardo linear z=0")

comparison("nonlin_ps_camb_z0.0.dat","r2xi_nonlin_camb_z0.0.dat","camb nonlinear z=0")

comparison("ps_z1.0.dat","r2xilin_z1.0.dat","Eduardo linear z=1")

comparison("nonlin_ps_z0.0.dat","r2xi_nonlin_z0.0.dat","Eduardo nonlinear z=0")

#For this one, the P and Xi are almost certainly not from the same model
#There are systematic shifts that look like Halofit_takahashi additions
comparison("nonlin_ps_z1.0.dat","r2xi_nonlin_z1.0.dat","Eduardo nonlinear z=1")
