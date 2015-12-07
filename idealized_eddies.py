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
from datetime import datetime, timedelta
from glob import glob
from tracpy.tracpy_class import Tracpy


def init(name, grid_filename, currents_filename):
    '''
    Initialize parameters for idealized eddy simulations. 
    Initialization for seeding drifters at all shelf model 
    grid points to be run forward.
    '''

    time_units = 'seconds since 0001-01-01  00:00:00'

    # horizontal_diffusivity project showed that relative dispersion did not
    # change between nsteps=25 and 50, but does between nsteps=5 and 25, and
    # interim numbers have not been tested yet.
    nsteps = 19  # to approximate the output timing of the TXLA model 25 

    # Number of steps to divide model output for outputting drifter location
    N = 4  # to approximate the output timing of the TXLA model # 5

    # Number of days
    ndays = 10

    # This is a forward-moving simulation
    ff = 1

    # Time between outputs
    tseas = 10800.0  # time between output in seconds
    ah = 0.
    av = 0.  # m^2/s

    # surface drifters
    z0 = 's'
    zpar = 29  # 30 layers

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # for periodic boundary conditions in the x direction
    doperiodic = 1

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    dostream = 0

    # Initialize Tracpy class
    tp = Tracpy(currents_filename, grid_filename, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps,
                N=N, ff=ff, ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar, time_units=time_units,
                usespherical=False, savell=False, doperiodic=doperiodic)

    # force grid reading
    tp._readgrid()

    # Start uniform array of drifters across domain using x,y coords
    x0 = tp.grid['xr'][1:-1,1:-1]
    y0 = tp.grid['yr'][1:-1,1:-1]

    return tp, x0, y0


def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('figures'):
        os.makedirs('figures')

    # Loop through all the simulations
    basepath = '/merrimack/raid/rob/Projects/shelfstrat/simulations/'
    # runs = ['shelfstrat_M2_1.49e-07_N2_1.00e-04_f_3.33e-05', 'shelfstrat_M2_1.92e-07_N2_1.00e-04_f_3.33e-05',
    #         'shelfstrat_M2_2.24e-07_N2_1.00e-04_f_5.00e-05', 'shelfstrat_M2_2.36e-07_N2_1.00e-04_f_3.33e-05',
    #         'shelfstrat_M2_2.89e-07_N2_1.00e-04_f_5.00e-05', 'shelfstrat_M2_3.33e-07_N2_1.00e-04_f_3.33e-05',
    #         'shelfstrat_M2_3.54e-07_N2_1.00e-04_f_5.00e-05', 'shelfstrat_M2_4.47e-07_N2_1.00e-04_f_1.00e-04',
    #         'shelfstrat_M2_5.00e-07_N2_1.00e-04_f_5.00e-05', 'shelfstrat_M2_5.77e-07_N2_1.00e-04_f_1.00e-04',
    #         'shelfstrat_M2_7.07e-07_N2_1.00e-04_f_1.00e-04', 'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04']

    # T update. For T=2.
    runs = ['shelfstrat_M2_4.47e-07_N2_1.00e-04_f_1.00e-04', 'shelfstrat_M2_5.00e-07_N2_1.00e-04_f_5.00e-05',
            'shelfstrat_M2_5.77e-07_N2_1.00e-04_f_1.00e-04', 'shelfstrat_M2_7.07e-07_N2_1.00e-04_f_1.00e-04',
            'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04']

    # list of start day for each simulation, in order. Calculated in Evernote.
    # startday = [7, 5, 4, 4, 3, 3, 2, 1, 2, 1, 1, 1]  # For T=3
    # startday = [10, 7, 5, 5, 3, 3, 3, 2, 2, 1, 1, 1]  # For T=4
    # startday = [23, 16, 11, 12, 8, 8, 6, 4, 4, 3, 2, 2]  # For T=10
    # startday = [46, 32, 21, 24, 15, 16, 12, 7, 8, 5, 4, 3]  # For T=20
    # startday = [12, 4]  # T=30
    
    # For updated T=2:
    startday = [16, 20, 12, 10, 7]

    for i, run in enumerate(runs):
        print run

        runname = run

        # Make symbolic links for run files from simulation in the main directory
        hisfileloc = basepath + run + '/shelfstrat_his.nc'
        grdfileloc = basepath + run + '/shelfstrat_grd.nc'

        # just do one simulation now
        f = netCDF.Dataset(hisfileloc)
        # t = f.variables['ocean_time']
        startdate = datetime(0001, 1, 1, 0, 0) + timedelta(days=startday[i])
        name = runname + '-startday' + str(startday[i])

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc') and \
           not os.path.exists('tracks/' + name + 'gc.nc'):

            # Read in simulation initialization
            tp, x0, y0 = init(name, grdfileloc, hisfileloc)

            xp, yp, zp, t, T0, U, V = tracpy.run.run(tp, startdate, x0, y0)

        # # Increment by 24 hours for next loop, to move through more quickly
        # date = date + timedelta(hours=24)



if __name__ == "__main__":
    run()    

