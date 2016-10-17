import sys
sys.path.insert(0,"src/")
import fastcorr
import numpy as np

#Load in a power spectrum
knl = np.genfromtxt("test_data/knl.txt")
pnl = np.genfromtxt("test_data/pnl.txt")

#Define the radial locations
NR = 1000
R = np.linspace(.1,120,NR)

#Call fastcorr
xi_nl = fastcorr.calc_corr(R,knl,pnl)

print xi_nl.shape
print xi_nl[:5]
