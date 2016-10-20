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
xi_nl = fastcorr.calc_corr(R,knl,pnl)

print xi_nl[:10]

import matplotlib.pyplot as plt
plt.loglog(R,R*R*xi_nl)
plt.xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
plt.ylabel(r"$R^2\xi(R)\ [{\rm Mpc^2}/h^2]$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
plt.show()

