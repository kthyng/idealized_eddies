'''
Functions to initialize various numerical experiments.

Make a new init_* for your application.

loc     Path to directory of grid and output files
nsteps  Number of steps to do between model outputs (iter in tracmass)
ndays   number of days to track the particles from start date
ff      ff=1 to go forward in time and ff=-1 for backward in time
date    Start date in datetime object
tseas   Time between outputs in seconds
ah      Horizontal diffusion in m^2/s. 
        See project values of 350, 100, 0, 2000. For -turb,-diffusion
av      Vertical diffusion in m^2/s.
do3d    for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
doturb  turbulence/diffusion flag. 
        doturb=0 means no turb/diffusion,
        doturb=1 means adding parameterized turbulence
        doturb=2 means adding diffusion on a circle
        doturb=3 means adding diffusion on an ellipse (anisodiffusion)
lon0    Drifter starting locations in x/zonal direction.
lat0    Drifter starting locations in y/meridional direction.
z0/zpar Then z0 should be an array of initial drifter depths. 
        The array should be the same size as lon0 and be negative
        for under water. Currently drifter depths need to be above 
        the seabed for every x,y particle location for the script to run.
        To do 3D but start at surface, use z0=zeros(ia.shape) and have
         either zpar='fromMSL'
        choose fromMSL to have z0 starting depths be for that depth below the base 
        time-independent sea level (or mean sea level).
        choose 'fromZeta' to have z0 starting depths be for that depth below the
        time-dependent sea surface. Haven't quite finished the 'fromZeta' case.
        Then: 
        set z0 to 's' for 2D along a terrain-following slice
         and zpar to be the index of s level you want to use (0 to km-1)
        set z0 to 'rho' for 2D along a density surface
         and zpar to be the density value you want to use
         Can do the same thing with salinity ('salt') or temperature ('temp')
         The model output doesn't currently have density though.
        set z0 to 'z' for 2D along a depth slice
         and zpar to be the constant (negative) depth value you want to use
        To simulate drifters at the surface, set z0 to 's' 
         and zpar = grid['km']-1 to put them in the upper s level
         z0='s' is currently not working correctly!!!
         In the meantime, do surface using the 3d set up option but with 2d flag set
xp      x-locations in x,y coordinates for drifters
yp      y-locations in x,y coordinates for drifters
zp      z-locations (depths from mean sea level) for drifters
t       time for drifter tracks
name    Name of simulation to be used for netcdf file containing final tracks

'''

import numpy as np
import os
import netCDF4 as netCDF
import pdb
import glob
from datetime import datetime, timedelta
from matplotlib.mlab import *
import tracpy
import time
from matplotlib import delaunay
import octant
from tracpy.tracpy_class import Tracpy

time_units = 'seconds since 0001-01-01'

def init(name, ndays, grid_filename, currents_filename):
    '''
    Initialization for seeding drifters at all shelf model grid points to be run
    forward.
    '''

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

    # Initial lon/lat locations for drifters
    # The following few lines aren't necessary because the grid cells are uniformly 1km res
    # # Start uniform array of drifters across domain using x,y coords
    # dx = 1000 # initial separation distance of drifters, in meters
    # xcrnrs = np.array([grid['xr'][1:-1,:].min(), grid['xr'][1:-1,:].max()])
    # ycrnrs = np.array([grid['yr'][1:-1,:].min(), grid['yr'][1:-1,:].max()])
    # X, Y = np.meshgrid(np.arange(xcrnrs[0], xcrnrs[1], dx), np.arange(ycrnrs[0], ycrnrs[1], dx))
    # lon0, lat0 = grid['basemap'](X, Y, inverse=True)

    # surface drifters
    z0 = 's'  
    zpar = 29 # 30 layers
    # zpar = grid['km']-1

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

    # return nsteps, N, ff, tseas, ah, av, x0, y0, \
    #         z0, zpar, do3d, doturb, grid, dostream, doperiodic
