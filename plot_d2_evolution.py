"""
Plot evolution of dispersion with increasing drifter start time.
"""

import numpy as np
from glob import glob
from matplotlib import colors
import cmocean
import matplotlib.pyplot as plt


# Read in dispersion calculations
locD2 = '/Volumes/COAWST/kristen/projects/idealized_eddies/calcs/D2/'
# locD2 = '/home/kthyng/projects/idealized_eddies/calcs/D2/'
Files = sorted(glob(locD2 + 'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04-startday*.npz'))

norm = colors.Normalize(0, 1)

# Plot
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
for i, File in enumerate(Files):
    d = np.load(File)
    D2 = d['D2']; days = d['days']; nnans = d['nnans']
    d.close()

    frac = i/20.
    color = cmocean.cm.cdom(norm(frac))

    ax.semilogy(days + i, D2, '--', color=color, lw=4, mec=None, alpha=0.6)


# shelfloc = '/rho/raid/home/kthyng/projects/shelf_transport/'
# fnameS = shelfloc + 'calcs/dispersion/hist/D2/r1/D2aveS.npz'
# D2aveS = np.load(fnameS)['D2aveS']
# dtracks = netCDF.Dataset(shelfloc + 'tracks/2004-01-01T00gc.nc')
# tp = dtracks.variables['tp']
# days2 = (tp-tp[0])/(3600.*24)
# dtracks.close()
# i25 = find(days2<=25)[-1]
# plt.semilogy(days2[:i25], D2aveS[70,45,:i25], '-', color='0.1', lw=4) # summer

# ax.set_xlim(0, 7)
