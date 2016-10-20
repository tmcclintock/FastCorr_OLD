import fastcorr
import time
import numpy as np

#Load in a power spectrum
knl = np.genfromtxt("../test_data/knl.txt")
pnl = np.genfromtxt("../test_data/pnl.txt")

#Define the radial locations
NR = 1000
R = np.linspace(0.1,200,NR)

#Call fastcorr
N,h = 1500,0.001

start = time.time()
xi_nl = fastcorr.calc_corr(R,knl,pnl,N,h)
end = time.time()
print "fastcorr xi0 time: %f"%(end-start)

start = time.time()
xi2_nl = fastcorr.calc_xi2(R,knl,pnl,N,h)
end = time.time()
print "fastcorr xi2 time: %f"%(end-start)

start = time.time()
xi4_nl = fastcorr.calc_xi4(R,knl,pnl,N,h)
end = time.time()
print "fastcorr xi2 time: %f"%(end-start)

#Properly scale the different parts
b,f = 2.0,0.75
xi_nl *= b**2 + 2./3.*b*f + 0.2*f**2
xi2_nl *= 4./3.*b*f + 4./7.*f**2
xi4_nl *= 8./35.*f**2

import matplotlib.pyplot as plt
plt.loglog(R,R*R*xi_nl,label=r"$R^2\xi_0$")
plt.loglog(R,-R*R*xi2_nl,label=r"$-R^2\xi_2$")
plt.loglog(R,R*R*xi4_nl,label=r"$R^2\xi_4$")
plt.ylim(1,1000)
plt.xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
plt.ylabel(r"$R^2\xi(R)\ [{\rm Mpc^2}/h^2]$",fontsize=24)
plt.legend().get_frame().set_alpha(0.5)
plt.subplots_adjust(bottom=0.15)
plt.show()

