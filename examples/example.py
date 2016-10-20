import fastcorr
import numpy as np

#Load in a power spectrum
knl = np.genfromtxt("../test_data/knl.txt")
pnl = np.genfromtxt("../test_data/pnl.txt")

#Define the radial locations
NR = 1000
R = np.linspace(10,200,NR)

b = 2
f = 0.75

#Call fastcorr
xi_nl = fastcorr.calc_corr(R,knl,pnl)*(b**2+2./3.*b*f+0.2*f**2)
#N,h = 3000,0.0002
xi2_nl = fastcorr.calc_xi2(R,knl,pnl)*(4./3.*b*f+4./7.*f**2)#,N,h)
xi4_nl = fastcorr.calc_xi4(R,knl,pnl)*(8./35.*f**2)#,N,h)
print xi4_nl[:10]

import matplotlib.pyplot as plt
plt.loglog(R,R*R*xi_nl)
plt.loglog(R,R*R*xi2_nl)
plt.loglog(R,-R*R*xi4_nl)

plt.show()

