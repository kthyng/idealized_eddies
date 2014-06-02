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

units = 'seconds since 0001-01-01'

def init():
    '''
    Initialization for seeding drifters at all shelf model grid points to be run
    forward.

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0). If date not input,
                name will be 'temp' 
        loc     Location of model output
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    # # Need to make a fake grid
    # loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    # grid = tracpy.inout.readgrid(loc)
    # # Replace entries
    currents_filename = 'ocean_his_0001.nc'
    grid_filename = 'grid.nc'
    # grid = tracpy.inout.readgrid(loc)
    # g = netCDF.Dataset('grid.nc')
    # o = netCDF.Dataset('ocean_his_0001.nc')
    # grid['xu'] = np.asfortranarray(g.variables['x_u'][:].T)
    # grid['xv'] = np.asfortranarray(g.variables['x_v'][:].T)
    # grid['yu'] = np.asfortranarray(g.variables['y_u'][:].T)
    # grid['yv'] = np.asfortranarray(g.variables['y_v'][:].T)
    # grid['xr'] = np.asfortranarray(g.variables['x_rho'][:].T)
    # grid['yr'] = np.asfortranarray(g.variables['y_rho'][:].T)
    # grid['xpsi'] = np.asfortranarray(g.variables['x_psi'][:].T)
    # grid['ypsi'] = np.asfortranarray(g.variables['y_psi'][:].T)

    # grid['lonr'], grid['latr'] = grid['basemap'](grid['xr'], grid['yr'], inverse=True)
    # grid['lonu'], grid['latu'] = grid['basemap'](grid['xu'], grid['yu'], inverse=True)
    # grid['lonv'], grid['latv'] = grid['basemap'](grid['xv'], grid['yv'], inverse=True)
    # grid['lonpsi'], grid['latpsi'] = grid['basemap'](grid['xpsi'], grid['ypsi'], inverse=True)
    # grid['lonr'] = np.asfortranarray(grid['lonr'])
    # grid['latr'] = np.asfortranarray(grid['latr'])
    # grid['lonu'] = np.asfortranarray(grid['lonu'])
    # grid['latu'] = np.asfortranarray(grid['latu'])
    # grid['lonv'] = np.asfortranarray(grid['lonv'])
    # grid['latv'] = np.asfortranarray(grid['latv'])
    # grid['lonpsi'] = np.asfortranarray(grid['lonpsi'])
    # grid['latpsi'] = np.asfortranarray(grid['latpsi'])

    # grid['mask'] = np.asfortranarray(g.variables['mask_rho'][:].T)
    # grid['pm'] = np.asfortranarray(g.variables['pm'][:].T)
    # grid['pn'] = np.asfortranarray(g.variables['pn'][:].T)

    # dxv = 1/grid['pm'] #.copy() # pm is 1/\Delta x at cell centers
    # dyu = 1/grid['pn'] #.copy() # pn is 1/\Delta y at cell centers

    # grid['dxdy'] = dyu*dxv
    # # # Already transposed
    # # grid['dxv'] = 0.5*(dxv[:-1,:]+dxv[1:,:])
    # # grid['dyu'] = 0.5*(dyu[:,:-1]+dyu[:,1:])
    # grid['dxv'] = 0.5*(dxv[:,:-1]+dxv[:,1:])
    # grid['dyu'] = 0.5*(dyu[:-1,:]+dyu[1:,:])

    # grid['imt'] = grid['xr'].shape[0]
    # grid['jmt'] = grid['xr'].shape[1]

    # grid['h'] = np.asfortranarray(g.variables['h'][:].T)
    # grid['sc_r'] = o.variables['s_w'][:] # sigma coords, 31 layers
    # grid['Cs_r'] = o.variables['Cs_w'][:] # stretching curve in sigma coords, 31 layers
    # grid['hc'] = o.variables['hc'][:]
    # grid['theta_s'] = o.variables['theta_s'][:]
    # grid['theta_b'] = o.variables['theta_b'][:]
    # grid['Vtransform'] = o.variables['Vtransform'][0]
    # grid['Vstretching'] = o.variables['Vstretching'][0]
    # grid['km'] = grid['sc_r'].shape[0]-1

    # mask2 = grid['mask'].copy()
    # kmt = np.ones((grid['imt'],grid['jmt']),order='f')*grid['km']
    # ind = (mask2==1)
    # ind[0:grid['imt']-1,:] = ind[1:grid['imt'],:]
    # mask2[ind] = 1
    # ind = (mask2==1)
    # ind[:,0:grid['jmt']-1] = ind[:,1:grid['jmt']]
    # mask2[ind] = 1
    # ind = (mask2==0)
    # kmt[ind] = 0
    # grid['kmt'] = kmt
    # grid['Y'], grid['X'] = np.meshgrid(np.arange(grid['jmt']),np.arange(grid['imt'])) # grid in index coordinates, without ghost cells
    # # Triangulation for grid space to curvilinear space
    # grid['tri'] = delaunay.Triangulation(grid['X'].flatten(),grid['Y'].flatten())
    # # Triangulation for curvilinear space to grid space
    # grid['trir'] = delaunay.Triangulation(grid['xr'].flatten(),grid['yr'].flatten())
    # grid['trirllrho'] = delaunay.Triangulation(grid['lonr'].flatten(),grid['latr'].flatten())
    # grid['zwt0'] = octant.depths.get_zw(grid['Vtransform'], grid['Vstretching'], grid['km']+1, 
    #                 grid['theta_s'], grid['theta_b'], 
    #                 grid['h'].T.copy(order='c'), 
    #                 grid['hc'], zeta=0, Hscale=3)
    # grid['zrt0'] = octant.depths.get_zrho(grid['Vtransform'], grid['Vstretching'], grid['km'], 
    #                 grid['theta_s'], grid['theta_b'], 
    #                 grid['h'].T.copy(order='c'), 
    #                 grid['hc'], zeta=0, Hscale=3)
    # # Change dzt to tracmass/fortran ordering
    # grid['zwt0'] = grid['zwt0'].T.copy(order='f')
    # grid['zrt0'] = grid['zrt0'].T.copy(order='f')
    # # this should be the base grid layer thickness that doesn't change in time because it 
    # # is for the reference vertical level
    # grid['dzt0'] = grid['zwt0'][:,:,1:] - grid['zwt0'][:,:,:-1]

    # g.close()

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
    zpar = grid['km']-1

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
                dtFromTracmass=dtFromTracmass, usespherical=False)

    # Start uniform array of drifters across domain using x,y coords
    x0 = tp.grid['xr'][1:-1,1:-1]
    y0 = tp.grid['yr'][1:-1,1:-1]

    return tp, x0, y0

    # return nsteps, N, ff, tseas, ah, av, x0, y0, \
    #         z0, zpar, do3d, doturb, grid, dostream, doperiodic
