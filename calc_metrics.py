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

T = 'other'
dofig = False

if T == 30:
    runs = ['shelfstrat_M2_5.00e-07_N2_1.00e-04_f_5.00e-05','shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04']
elif T == 40:
    # just the base case
    runs = ['shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04','shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04',
            'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04','shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04',
            'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04']
elif T == 2:
    # T update. For T=2.
    runs = ['shelfstrat_M2_4.47e-07_N2_1.00e-04_f_1.00e-04', 'shelfstrat_M2_5.00e-07_N2_1.00e-04_f_5.00e-05',
            'shelfstrat_M2_5.77e-07_N2_1.00e-04_f_1.00e-04', 'shelfstrat_M2_7.07e-07_N2_1.00e-04_f_1.00e-04',
            'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04']
elif T == 'other':
    startday = np.arange(1, 24)
    runname = 'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04'
    runs = [runname for i in range(len(startday))]
else:
    runs = ['shelfstrat_M2_1.49e-07_N2_1.00e-04_f_3.33e-05', 'shelfstrat_M2_1.92e-07_N2_1.00e-04_f_3.33e-05',
            'shelfstrat_M2_2.24e-07_N2_1.00e-04_f_5.00e-05', 'shelfstrat_M2_2.36e-07_N2_1.00e-04_f_3.33e-05',
            'shelfstrat_M2_2.89e-07_N2_1.00e-04_f_5.00e-05', 'shelfstrat_M2_3.33e-07_N2_1.00e-04_f_3.33e-05',
            'shelfstrat_M2_3.54e-07_N2_1.00e-04_f_5.00e-05', 'shelfstrat_M2_4.47e-07_N2_1.00e-04_f_1.00e-04',
            'shelfstrat_M2_5.00e-07_N2_1.00e-04_f_5.00e-05', 'shelfstrat_M2_5.77e-07_N2_1.00e-04_f_1.00e-04',
            'shelfstrat_M2_7.07e-07_N2_1.00e-04_f_1.00e-04', 'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04']

if not T == 'other':
    # list of start day for each simulation, in order. Calculated in Evernote.
    if T == 3:
        startday = [7, 5, 4, 4, 3, 3, 2, 1, 2, 1, 1, 1]  # For T=3
    elif T == 4:
        startday = [10, 7, 5, 5, 3, 3, 3, 2, 2, 1, 1, 1]  # For T=4
    elif T == 10:
        startday = [23, 16, 11, 12, 8, 8, 6, 4, 4, 3, 2, 2]  # For T=10
    elif T == 20:
        startday = [46, 32, 21, 24, 15, 16, 12, 7, 8, 5, 4, 3]  # For T=20
    elif T == 30:
        startday = [12, 4]
    elif T == 40:
        startday = [3, 4, 5, 6, 7]  # only base case, more starting days
    elif T == 2:
        startday = [16, 20, 12, 10, 7]

grdfileloc = basepath + runs[0] + '/shelfstrat_grd.nc'
# grid = tracpy.inout.readgrid(grdfileloc, usespherical=False)

# Set up figure
if dofig:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)

for i, run in enumerate(runs):

    fname = 'calcs/D2/' + run + '-startday' + str(startday[i]).zfill(2) + '.npz'
    print fname

    if os.path.exists(fname):
        d = np.load(fname)
        D2 = d['D2']; days = d['days']; nnans = d['nnans']
        d.close()

    else:

        name = run + '-startday' + str(startday[i])

        # Use drifters starting once eddies are developed in all simulations: 20 days
        d = netCDF.Dataset(glob('tracks/' + name + '*.nc')[0])

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

        if dofig:
            # plot
            if '3.33e-05' in run: # row 3
                color = 'g' # '0.4'
            elif '5.00e-05' in run: # row 2
                color = 'darkorange' # '0.1'
            elif '1.00e-04' in run: # row 1
                color = 'r'
            # print run
            # dalpha = (1-.25)/4.
            # alpha = 0.25+dalpha*np.mod(i,5)
            linestyles = ['-', '--', '-.', ':']#, 'o:']
            if ('M2_3.33e-07' in run) or ('M2_5.00e-07' in run) or ('M2_1.00e-06' in run):
                linestyle = linestyles[0]
            elif ('M2_2.36e-07' in run) or ('M2_3.54e-07' in run) or ('M2_7.07e-07' in run):
                linestyle = linestyles[1]
            elif ('M2_1.92e-07' in run) or ('M2_2.89e-07' in run) or ('M2_5.77e-07' in run):
                linestyle = linestyles[2]
            elif ('M2_1.49e-07' in run) or ('M2_2.24e-07' in run) or ('M2_4.47e-07' in run):
                linestyle = linestyles[3]
            ax.semilogy(days, D2, linestyle, color=color, lw=4, mec=None, alpha=0.6)

    # for ipair in pairs:
    #     dist = np.sqrt((xg[ipair[0],:] - xg[ipair[1],:])**2 + (yg[ipair[0],:] - yg[ipair[1],:])**2)
    #     if np.sum(abs(np.diff(dist))>10)>1:
    #         print ipair
    #     if np.sum(dist>200)>1:
    #         print ipair

if dofig:
    # ax.set_xlabel('Days 20 to 30')
    ax.set_ylabel('Mean squared separation distance [km$^2\!$]')
    ax.set_xlim(0, 10)

    ax.text(0.05, 0.95, 'Row 1', color='r', transform=ax.transAxes)
    ax.text(0.05, 0.92, 'Row 2', color='orange', transform=ax.transAxes)
    ax.text(0.05, 0.89, 'Row 3', color='g', transform=ax.transAxes)

    fig.savefig('figures/D2-T' + str(T) + '.png', bbox_inches='tight')
    fig.savefig('figures/D2-T' + str(T) + '.pdf', bbox_inches='tight')

    # add on TXLA shelf eddy region comparison
    shelfloc = '/rho/raid/home/kthyng/projects/shelf_transport/'
    fnameS = shelfloc + 'calcs/dispersion/hist/D2/r1/D2aveS.npz'
    D2aveS = np.load(fnameS)['D2aveS']
    dtracks = netCDF.Dataset(shelfloc + 'tracks/2004-01-01T00gc.nc')
    tp = dtracks.variables['tp']
    days2 = (tp-tp[0])/(3600.*24)
    dtracks.close()
    i25 = find(days2<=25)[-1]
    plt.semilogy(days2[:i25], D2aveS[70,45,:i25], '-', color='0.3', lw=4) # summer

    fig.savefig('figures/D2-withTXLA-T' + str(T) + '.png', bbox_inches='tight')
    fig.savefig('figures/D2-withTXLA-T' + str(T) + '.pdf', bbox_inches='tight')

# if __name__ == "__main__":
#     run()    
