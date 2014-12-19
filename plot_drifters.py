'''
Plot a set of drifters in map view, in time.
'''

import matplotlib.pyplot as plt
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
from matplotlib.mlab import find
import pdb
import numpy as np
import matplotlib as mpl
import os
import tracpy.calcs
import glob

mpl.rcParams.update({'font.size': 20})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

# whether to do tails on drifters or not (don't with low decimation)
dotails = False # True or False
donewtails = True # for when not doing old tails

# Read in drifter tracks
dd = 1 # 500 # drifter decimation
startdate = '0001-01-21T00' #14T00'
# runs = ['tracks/shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04']
# runs = glob.glob('tracks/shelfstrat_M2_3.33e-07_N2_1.00e-04_f_7.45*')
# runs = glob.glob('tracks/shelfstrat_M2_1.00e-06*')
# runs row 1
base = 'tracks/shelfstrat_M2_2.00e-07_N2_1.00e-04_f_'
runs = [base + '2.00e-05', base + '2.83e-05', base + '3.46e-05', base + '4.47e-05', base + '6.32e-05']
# # runs row 2
# base = 'tracks/shelfstrat_M2_3.33e-07_N2_1.00e-04_f_'
# runs = [base + '3.33e-05', base + '4.71e-05', base + '5.77e-05', base + '7.45e-05', base + '1.05e-04']
# # runs row 3
# base = 'tracks/shelfstrat_M2_5.00e-07_N2_1.00e-04_f_'
# runs = [base + '5.00e-05', base + '7.07e-05', base + '8.66e-05', base + '1.12e-04', base + '1.58e-04']
# # runs row 4
# base = 'tracks/shelfstrat_M2_1.00e-06_N2_1.00e-04_f_'
# runs = [base + '1.00e-04', base + '1.41e-04', base + '1.73e-04', base + '2.24e-04', base + '3.16e-04']

for run in runs:
    # pdb.set_trace()

    # read in grid
    basepath = '/merrimack/raid/rob/Projects/shelfstrat/simulations/'
    grdfileloc = basepath + run.split('/')[-1] + '/shelfstrat_grd.nc'
    hisfileloc = basepath + run.split('/')[-1] + '/shelfstrat_his.nc'
    grid = tracpy.inout.readgrid(grdfileloc, usebasemap=False, usespherical=False)

    d = netCDF.Dataset(run + '/' + startdate + 'gc.nc')
    xg = d.variables['xg'][::dd,:]
    yg = d.variables['yg'][::dd,:]


    ind = abs(np.diff(xg))>3
    for j in xrange(ind.shape[0]):
        icrits = find(ind[j,:])
        xg[j,icrits] = np.nan
    # # Fix x position of drifters for periodic boundaries in the x direction
    # xg = np.unwrap((xg*(2*np.pi)/256.), axis=-1)*(256./(2*np.pi))

    # # Need to first eliminate the interpolated points between ends in the periodic direction
    # # Identify these points by having large jumps for two consecutive x-direction positions
    # ind = abs(np.diff(xg, axis=1))>3
    # for j in xrange(ind.shape[0]):
    #     doloop = True
    #     icrits = (find(ind[j,1:] & ind[j,:-1]) + 1).astype(list) # index of weird interpolated point
    #     if icrits.size == 0: # if this exists for this drifter
    #         doloop = False
    #     while doloop:
    #         icrit = icrits[0]
    #         xg[j,icrit] = .5*(xg[j,icrit-1] + xg[j,icrit+1])
    #         # update indices in case any of these were in a row
    #         ind[j,:] = abs(np.diff(xg[j,:]))>3
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

    # pdb.set_trace()


    xp = xg.copy()*1000; yp = yg.copy()*1000
    # ind = (xg == -1)
    # xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
    # xp[ind] = np.nan; yp[ind] = np.nan
    tp = d.variables['tp'][0,:]
    # d.close()

    # Change to km from m
    xp /= 1000.
    yp /= 1000.

    # # txla output
    # nc = netCDF.Dataset(loc)
    # datestxla = netCDF.num2date(nc.variables['ocean_time'][:], nc.variables['ocean_time'].units)

    if not dotails:
        
        # Find indices of drifters by starting depth
        depthp = tracpy.calcs.Var(xg[:,0], yg[:,0], tp, 'h', netCDF.Dataset(hisfileloc)) # starting depths of drifters
        # near-shore
        # ind10 = depthp<=10
        ind20 = depthp<=20
        ind50 = (depthp>20)*(depthp<=50)
        ind100 = (depthp>50)*(depthp<=100)
        # offshore
        ind500 = (depthp>100)*(depthp<=500)
        ind3500 = depthp>500

        # colors for drifters
        rgb = plt.cm.get_cmap('winter_r')(np.linspace(0,1,6))[:-1,:3] # skip last entry where it levels off in lightness

        # to plot colorbar
        gradient = np.linspace(0, 1, 6)[:-1]
        gradient = np.vstack((gradient, gradient))

        ms = 4 #1.5 # markersize


    # Plot drifters, starting 5 days into simulation
    # 2 days for one tail, 3 days for other tail
    # t = tp-tp[0]
    days = (tp-tp[0])/(3600.*24)
    dates = netCDF.num2date(tp, d.variables['tp'].units)
    # pdb.set_trace()
    # Find indices relative to present time
    i5daysago = 0 # keeps track of index 5 days ago
    i2daysago = find(days>=3)[0] # index for 2 days ago, which starts as 3 days in
    if dotails:
        i5days = find(days>=5)[0] # index for 5 days in
    else:
        i5days = 0 # start at the beginning
    nt = tp.size # total number of time indices
    # for i in np.arange(0,nt+1,5):

    dirname = 'figures/drifters/dd' + str(dd) + '/' + run.split('/')[-1] + '/' + startdate
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    di = 4
    for i in np.arange(i5days,nt,di): # only plot for every circulation model output
        # pdb.set_trace()
    # for i in np.arange(i5days,nt+1,5):
        # if i==545:
        #     pdb.set_trace()

        fname = dirname + '/' + str(dates[i])[3:-6].replace(' ', '-') + '.png'
        # if not dotails:
        #     itxla = np.where(datestxla==dates[i])[0][0] # find time index to use for model output
        #     salt = nc.variables['salt'][itxla,-1,:,:] # surface salinity

        if os.path.exists(fname):
            # Update indices
            i5daysago += di
            i2daysago += di
            continue

        # Plot background
        fig = plt.figure(figsize=(14,5.5))
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.07, bottom=0.12, right=0.99, top=1.0)

        if dotails:
            # Plot 5 days ago to 2 days ago
            ax.plot(xp[:,i5daysago:i2daysago].T, yp[:,i5daysago:i2daysago].T, color='0.6', lw=2)

            # Plot 0-2 day tail
            ax.plot(xp[:,i2daysago:i].T, yp[:,i2daysago:i].T, color='0.3', lw=3)

            # Plot drifter locations
            ax.plot(xp[:,i].T, yp[:,i].T, 'o', color='r', ms=10)

        else:

            if donewtails:

                istart = i-(24./3*4*0.25) # when to start plotting lines back in time
                if istart<0: istart=0

                if i==0: # plot dots in starting locations

                    # Plot drifter locations
                    ax.plot(xp[ind20,i].T, yp[ind20,i].T, '.', color=rgb[0,:], ms=0.6)
                    ax.plot(xp[ind50,i].T, yp[ind50,i].T, '.', color=rgb[1,:], ms=0.6)
                    ax.plot(xp[ind100,i].T, yp[ind100,i].T, '.', color=rgb[2,:], ms=0.6)
                    ax.plot(xp[ind500,i].T, yp[ind500,i].T, '.', color=rgb[3,:], ms=0.6)
                    ax.plot(xp[ind3500,i].T, yp[ind3500,i].T, '.', color=rgb[4,:], ms=0.6)

                else:

                    # Plot drifter locations
                    ax.plot(xp[ind20,istart:i+1].T, yp[ind20,istart:i+1].T, '-', color=rgb[0,:], lw=0.5)
                    ax.plot(xp[ind50,istart:i+1].T, yp[ind50,istart:i+1].T, '-', color=rgb[1,:], lw=0.5)
                    ax.plot(xp[ind100,istart:i+1].T, yp[ind100,istart:i+1].T, '-', color=rgb[2,:], lw=0.5)
                    ax.plot(xp[ind500,istart:i+1].T, yp[ind500,istart:i+1].T, '-', color=rgb[3,:], lw=0.5)
                    ax.plot(xp[ind3500,istart:i+1].T, yp[ind3500,istart:i+1].T, '-', color=rgb[4,:], lw=0.5)

            else:

                # Plot drifter locations
                # ax.plot(xp[:,i].T, yp[:,i].T, 'o', color='g', ms=2, mec='None')
                ax.plot(xp[ind20,i].T, yp[ind20,i].T, 'o', color=rgb[0,:], ms=ms, mec='None')
                ax.plot(xp[ind50,i].T, yp[ind50,i].T, 'o', color=rgb[1,:], ms=ms, mec='None')
                ax.plot(xp[ind100,i].T, yp[ind100,i].T, 'o', color=rgb[2,:], ms=ms, mec='None')
                ax.plot(xp[ind500,i].T, yp[ind500,i].T, 'o', color=rgb[3,:], ms=ms, mec='None')
                ax.plot(xp[ind3500,i].T, yp[ind3500,i].T, 'o', color=rgb[4,:], ms=ms, mec='None')

            # # Overlay surface salinity
            # ax.contour(grid['xr'].T, grid['yr'].T, salt, [33], colors='0.1', zorder=12, linewidths=2)

        # pdb.set_trace()
        # plt.axis('tight')
        ax.set_xlim(1,255)
        ax.set_ylim(1,90)#129)
        ax.set_aspect(1.0)
        ax.set_ylabel('Across channel [km]')
        ax.set_xlabel('Along channel [km]')

        # Time
        ax.text(0.0, -0.115, 'Days ' + str(days[i]), transform=ax.transAxes, fontsize=20)
        # ax.text(0.075, 0.95, dates[i].isoformat()[:-6], transform=ax.transAxes, fontsize=20)

        # Drifter legend
        if dotails:
            ax.plot(0.0895, 0.9, 'or', ms=10, transform=ax.transAxes) # drifter head
            ax.plot([0.075, 0.1], [0.875, 0.875], '0.3', lw=3, transform=ax.transAxes) # drifter tail #1
            ax.plot([0.075, 0.1], [0.85, 0.85], '0.5', lw=2, transform=ax.transAxes) # drifter tail #2
            ax.text(0.125, 0.89, 'Drifter location', color='r', transform=ax.transAxes, fontsize=16)
            ax.text(0.125, 0.866, '2 days prior', color='0.3', transform=ax.transAxes, fontsize=16)
            ax.text(0.125, 0.842, '5 days prior', color='0.5', transform=ax.transAxes, fontsize=16)
        else:
            cax = fig.add_axes([0.79, 0.04, 0.175, 0.02])
            cax.imshow(gradient, aspect='auto', interpolation='none', cmap=plt.get_cmap('winter_r'))
            cax.tick_params(axis='y', labelleft=False, left=False, right=False)
            cax.tick_params(axis='x', top=False, bottom=False, labelsize=15)
            cax.set_xticks(np.arange(-0.5, 5, 1.0))
            cax.set_xticklabels(('0', '20', '50', '100', '500', '3500'))
            cax.set_title('Initial drifter depth [m]', fontsize=14)
        # pdb.set_trace()

            # # legend for contour
            # ax.plot([0.075, 0.1], [0.81, 0.81], '0.1', lw=2, transform=ax.transAxes)
            # ax.text(0.123, 0.802, '33 salinity contour', color='0.1', transform=ax.transAxes, fontsize=16)


        # Update indices
        i5daysago += di
        i2daysago += di

        fig.savefig(fname)

        plt.close()