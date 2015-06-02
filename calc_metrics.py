'''
Calculate metrics of drifters starting from specific areas.
'''

import numpy as np
import pdb
from matplotlib.mlab import find
import netCDF4 as netCDF
from scipy import ndimage
import time
from glob import glob
import tracpy
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import tracpy.calcs

# mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 16})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


basepath = '/merrimack/raid/rob/Projects/shelfstrat/simulations/'

# get simulation names to loop through
sims = glob(basepath + 'shelfstrat*')

grdfileloc = basepath + sims[0] + '/shelfstrat_grd.nc'
# grid = tracpy.inout.readgrid(grdfileloc, usespherical=False)

# Set up figure
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)

for i, sim in enumerate(sims):

    fname = 'calcs/D2/' + sim.split('/')[-1] + '.npz'
    print fname

    if os.path.exists(fname):

        d = np.load(fname)
        D2 = d['D2']; days = d['days']; nnans = d['nnans']
        d.close()

    else:

        # Use drifters starting once eddies are developed in all simulations: 20 days
        d = netCDF.Dataset(glob('tracks/' + sim.split('/')[-1] + '/*')[0])

        tp = d.variables['tp'][0,:]
        days = (tp-tp[0])/(3600.*24)
        xg = d.variables['xg'][:]; yg = d.variables['yg'][:]

        # Only use drifters that start within 60 km of the shore in the y direction
        inds = yg[:,0]<=60 # 60km is grid index 60
        xg = xg[inds,:]; yg = yg[inds,:]

        # Unwrap locations in the periodic direction
        # Do this a first time to get all the easy unwrapping done so it won't show up
        # in the next logic tests
        xg = np.unwrap((xg*(2*np.pi)/256.), axis=-1)*(256./(2*np.pi))

        # # Need to first eliminate the interpolated points between ends in the periodic direction
        # # Identify these points by having large jumps for two consecutive x-direction positions
        # print 'preprocessing...'
        # ind = abs(np.diff(xg, axis=1))>5
        # for j in xrange(ind.shape[0]):
        #     doloop = True
        #     icrits = (find(ind[j,1:] & ind[j,:-1]) + 1).astype(list) # index of weird interpolated point
        #     if icrits.size == 0: # if this exists for this drifter
        #         doloop = False
        #     while doloop:
        #         if 142 in icrits:
        #             pdb.set_trace()
        #         icrit = icrits[0]
        #         print icrit
        #         xg[j,icrit] = .5*(xg[j,icrit-1] + xg[j,icrit+1])
        #         # update indices in case any of these were in a row
        #         ind[j,:] = abs(np.diff(xg[j,:]))>5
        #         icrits = find(ind[j,1:] & ind[j,:-1]) + 1 # index of weird interpolated point
        #         # test to see if this worked by seeing if the same critical index is still there
        #         if (icrits.size > 0) and (icrits[0] == icrit): # then it didn't work, try a different approach
        #             if xg[j,icrit+1]>250: # trying to wrap left to right
        #                 xg[j,icrit] = xg[j,icrit-1] - (256 - xg[j,icrit+1])
        #             elif xg[j,icrit+1]<100: # trying to wrap right to left
        #                 xg[j,icrit] = xg[j,icrit-1] + (xg[j,icrit+1])
        #             # unwrap again
        #             xg[j,:] = np.unwrap((xg[j,:]*(2*np.pi)/256.), axis=-1)*(256./(2*np.pi))
                
        #         if icrits.size == 0:
        #             doloop = False

        d.close()

        # Run calculation
        print 'relative dispersion calculation...'
        D2, nnans, pairs = tracpy.calcs.rel_dispersion(xg*1000, yg*1000, r=[0.95, 1.05], squared=True, spherical=False)

        # Save
        np.savez(fname, D2=D2, nnans=nnans, pairs=pairs, days=days)

    # plot
    if '2.00e-07' in sim: # row 4
        color = 'b' # '0.7'
    elif '3.33e-07' in sim: # row 3
        color = 'g' # '0.4'
    elif '5.00e-07' in sim: # row 2
        color = 'darkorange' # '0.1'
    elif '1.00e-06' in sim: # row 1
        color = 'r'
    # print sim
    # dalpha = (1-.25)/4.
    # alpha = 0.25+dalpha*np.mod(i,5)
    linestyles = ['-', '--', '-.', ':', 'o:']
    dds = [1,1,1,1,20]
    ax.semilogy(days[::dds[np.mod(i,5)]], D2[::dds[np.mod(i,5)]], linestyles[np.mod(i,5)], color=color, lw=4, mec=None, alpha=0.6)

    # for ipair in pairs:
    #     dist = np.sqrt((xg[ipair[0],:] - xg[ipair[1],:])**2 + (yg[ipair[0],:] - yg[ipair[1],:])**2)
    #     if np.sum(abs(np.diff(dist))>10)>1:
    #         print ipair
    #     if np.sum(dist>200)>1:
    #         print ipair

ax.set_xlabel('Days 20 to 30')
ax.set_ylabel('Mean squared separation distance [km$^2\!$]')

ax.text(0.05, 0.95, 'Row 1', color='r', transform=ax.transAxes)
ax.text(0.05, 0.92, 'Row 2', color='orange', transform=ax.transAxes)
ax.text(0.05, 0.89, 'Row 3', color='g', transform=ax.transAxes)
ax.text(0.05, 0.86, 'Row 4', color='b', transform=ax.transAxes)

fig.savefig('figures/D2.png', bbox_inches='tight')
fig.savefig('figures/D2.pdf', bbox_inches='tight')

# add on TXLA shelf eddy region comparison
shelfloc = '/rho/raid/home/kthyng/projects/shelf_transport/'
fnameS = shelfloc + 'calcs/dispersion/hist/D2/r1/D2aveS.npz'
D2aveS = np.load(fnameS)['D2aveS']
dtracks = netCDF.Dataset(shelfloc + 'tracks/2004-01-01T00gc.nc')
tp = dtracks.variables['tp']
days2 = (tp-tp[0])/(3600.*24)
dtracks.close()
i10 = find(days2<=10)[-1]
plt.semilogy(days2[:i10], D2aveS[70,45,:i10], '-', color='0.3', lw=4) # summer

fig.savefig('figures/D2-withTXLA.png', bbox_inches='tight')
fig.savefig('figures/D2-withTXLA.pdf', bbox_inches='tight')

# if __name__ == "__main__":
#     run()    
