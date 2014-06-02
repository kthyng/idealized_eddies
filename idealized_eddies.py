'''
Script to run drifters at 1km initial spacing daily forward for 30 days.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
from glob import glob

# npieces = 12 # number of pieces to divide starting locations for drifters into, in x direction

def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('figures'):
        os.makedirs('figures')

    # Loop through all the simulations
    # basepath = '/merrimack/raid/rob/Projects/shelfstrat/simulations-unbalanced/'
    # runs = glob(basepath + 'shelfstrat*')
    # # runs = [basepath + 'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04']
    basepath = '/merrimack/raid/rob/Projects/shelfstrat/simulations-balanced/'
    # runs = glob(basepath + 'shelfstrat*')
    runs = [basepath + 'shelfstrat_M2_8.94e-08_N2_1.00e-04_f_2.00e-05']

    for run in runs:

        runname = 'balanced/' + run.split('/')[-1]
        # runname = 'unbalanced/' + run.split('/')[-1]

        # make directory for output
        if not os.path.exists('tracks/' + runname):
            os.makedirs('tracks/' + runname)

        # Make symbolic links for run files from simulation in the main directory
        hisfileloc = run + '/shelfstrat_his.nc'
        grdfileloc = run + '/shelfstrat_grd.nc'
        hisfiledes = 'ocean_his_0001.nc' # 'tracks/' + runname + '/ocean_his_0001.nc'
        grdfiledes = 'grid.nc' # 'tracks/' + runname + '/grid.nc'
        hiscmd = 'ln -sf ' + hisfileloc + ' ' + hisfiledes
        grdcmd = 'ln -sf ' + grdfileloc + ' ' + grdfiledes
        # pdb.set_trace()
        os.system(hiscmd)
        os.system(grdcmd)

        # Start drifters daily in the simulations
        overallstartdate = datetime(0001, 1, 1, 0, 0)
        overallstopdate = datetime(0001, 1, 31, 0, 0)

        date = overallstartdate

        # Start from the beginning and add days on for loop
        # keep running until we hit the next month
        while date < overallstopdate:

            name = date.isoformat()[0:13]

            # If the particle trajectories have not been run, run them
            if not os.path.exists('tracks/' + runname + '/' + name + '.nc') and \
                not os.path.exists('tracks/' + runname + '/' + name + 'gc.nc'):

                # Read in simulation initialization
                ndays = (overallstopdate-date).days
                tp, x0, y0 = init.init()
                # nstep, N, ff, tseas, ah, av, lon0, lat0, z0, zpar, do3d, doturb, \
                #         grid, dostream, doperiodic = init.init()
                # pdb.set_trace()
                # Run tracpy
                # Save directly to grid coordinates
                # lonp, latp, zp, t, grid \
                #     = tracpy.run.run(['ocean_his_0001.nc',''], 
                #                         nstep, ndays, ff, date, tseas, ah, av, \
                #                         lon0, lat0, z0, zpar, do3d, doturb, runname + '/' + name, N=N,  \
                #                         grid=grid, dostream=dostream, doperiodic=doperiodic, savell=False)
                xp, yp, zp, t, T0, U, V = tracpy.run.run(tp, date, x0, y0)

            # # If basic figures don't exist, make them
            # if not os.path.exists('figures/' + name + '*.png'):

                # # Read in and plot tracks
                # d = netCDF.Dataset('tracks/' + name + '.nc')
                # lonp = d.variables['lonp'][:]
                # latp = d.variables['latp'][:]
                # # tracpy.plotting.tracks(lonp, latp, name, grid=grid)
                # # tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
                # d.close()
                # # # Do transport plot
                # tracpy.plotting.transport(name='all_f/N=5_dx=8/25days', fmod=date.isoformat()[0:13], 
                #     extraname=date.isoformat()[0:13], 
                #     Title='Transport on Shelf, for a week from ' + date.isoformat()[0:13], dmax=1.0)

            # Increment by 24 hours for next loop, to move through more quickly
            # nh = nh + 24
            date = date + timedelta(hours=24)

        # # Do transport plot
        # tracpy.plotting.transport(name='all_f/N=5_dx=8/25days', fmod=startdate.isoformat()[0:7] + '*', 
        #     extraname=startdate.isoformat()[0:7], Title='Transport on Shelf', dmax=1.0)


if __name__ == "__main__":
    run()    
