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
from tracpy.tracpy_class import Tracpy


def init(name, ndays, grid_filename, currents_filename):
    '''
    Initialize parameters for idealized eddy simulations. 
    Initialization for seeding drifters at all shelf model 
    grid points to be run forward.
    '''

    time_units = 'seconds since 0001-01-01'

    # horizontal_diffusivity project showed that relative dispersion did not
    # change between nsteps=25 and 50, but does between nsteps=5 and 25, and
    # interim numbers have not been tested yet.
    nsteps = 19 # to approximate the output timing of the TXLA model 25 

    # Number of steps to divide model output for outputting drifter location
    N = 4 # to approximate the output timing of the TXLA model # 5

    # Number of days
    # ndays = 30

    # This is a forward-moving simulation
    ff = 1 

    # Time between outputs
    tseas = 10800.0 # time between output in seconds
    ah = 0.
    av = 0. # m^2/s

    # surface drifters
    z0 = 's'  
    zpar = 29 # 30 layers

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
                usespherical=False, savell=False, doperiodic=1)

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
    runs = glob(basepath + 'shelfstrat*')

    for run in runs:

        runname = run.split('/')[-1]

        # make directory for output
        if not os.path.exists('tracks/' + runname):
            os.makedirs('tracks/' + runname)

        # pdb.set_trace()

        # Make symbolic links for run files from simulation in the main directory
        hisfileloc = run + '/shelfstrat_his.nc'
        grdfileloc = run + '/shelfstrat_grd.nc'
        # hisfiledes = 'ocean_his_0001.nc' # 'tracks/' + runname + '/ocean_his_0001.nc'
        # grdfiledes = 'grid.nc' # 'tracks/' + runname + '/grid.nc'
        # hiscmd = 'ln -sf ' + hisfileloc + ' ' + hisfiledes
        # grdcmd = 'ln -sf ' + grdfileloc + ' ' + grdfiledes
        # pdb.set_trace()
        # os.system(hiscmd)
        # os.system(grdcmd)

        # just do one simulation now
        f = netCDF.Dataset(hisfileloc)
        t = f.variables['ocean_time'][:]
        days = t[-1]/86400. # total number of days for simulation
        daysdate = datetime(0001, 1, 1, 0, 0) + timedelta(days=int(days)) # end date of ocean simulation
        date = daysdate - timedelta(days=7) # date to start drifter simulation

        # # Start drifters daily in the simulations
        # overallstartdate = datetime(0001, 1, 1, 0, 0)
        # overallstopdate = datetime(0001, 1, 31, 0, 0)

        # date = overallstartdate

        # # Start from the beginning and add days on for loop
        # # keep running until we hit the next month
        # while date < overallstopdate:

        name = date.isoformat()[0:13]

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + runname + '/' + name + '.nc') and \
            not os.path.exists('tracks/' + runname + '/' + name + 'gc.nc'):

            # Read in simulation initialization
            ndays = 7 #(overallstopdate-date).days
            tp, x0, y0 = init(runname + '/' + name, ndays, grdfileloc, hisfileloc)

            xp, yp, zp, t, T0, U, V = tracpy.run.run(tp, date, x0, y0)

        # # Increment by 24 hours for next loop, to move through more quickly
        # date = date + timedelta(hours=24)



if __name__ == "__main__":
    run()    
