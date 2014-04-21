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


def select_pts():
    '''
    Select out which initial drifter locations to use in calculation.
    '''

    # Distance in indices in y direction for each strip
    dy = 5
    # Number of grid indices to use, total, in y direction
    ny = 80
    nx = 258

    strips = [] # to store points

    for i in xrange(ny/dy):
        strips.append([np.arange(i*dy,i*dy+dy), np.arange(nx)])

    return strips


def select_times(dosim):
    '''
    return different time periods for calculation based on which simulation using.
    '''

    # stores the indices for the ramp up period in 0 and for steady timing in 1
    tinds = []
    itmax = 901 # max number of time indices

    if dosim == 'shelfstrat_M2_5.00e-07_N2_1.00e-04_f_5.00e-05':
        itramp = 183; itsteady = 330
    elif dosim == 'shelfstrat_M2_3.54e-07_N2_1.00e-04_f_5.00e-05':
        itramp = 640; itsteady = itmax
    elif dosim == 'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04':
        itramp = 112; itsteady = 218
    elif dosim == 'shelfstrat_M2_7.07e-07_N2_1.00e-04_f_1.00e-04':
        itramp = 160; itsteady = 454
    elif dosim == 'shelfstrat_M2_4.47e-07_N2_1.00e-04_f_1.00e-04':
        itramp = 256; itsteady = 551
    elif dosim == 'shelfstrat_M2_3.16e-07_N2_1.00e-04_f_1.00e-04':
        itramp = 362; itsteady = 768

    tinds.append(np.arange(itramp,itsteady)) # ramp up period
    tinds.append(np.arange(itsteady,itmax+1)) # ramp up period

    return tinds

# def select_times(dosim):
#     '''
#     return different time periods for calculation based on which simulation using.
#     '''

#     # stores the indices for the ramp up period in 0 and for steady timing in 1
#     tinds = []
#     itmax = 901 # max number of time indices

#     if dosim == 'shelfstrat_M2_5.00e-07_N2_1.00e-04_f_5.00e-05':
#         itramp = 183; itsteady = 330
#     elif dosim == 'shelfstrat_M2_3.54e-07_N2_1.00e-04_f_5.00e-05':
#         itramp = 640; itsteady = itmax
#     elif dosim == 'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04':
#         itramp = 112; itsteady = 218
#     elif dosim == 'shelfstrat_M2_7.07e-07_N2_1.00e-04_f_1.00e-04':
#         itramp = 160; itsteady = 454
#     elif dosim == 'shelfstrat_M2_4.47e-07_N2_1.00e-04_f_1.00e-04':
#         itramp = 256; itsteady = 551
#     elif dosim == 'shelfstrat_M2_3.16e-07_N2_1.00e-04_f_1.00e-04':
#         itramp = 362; itsteady = 768

#     tinds.append(np.arange(itramp,itsteady)) # ramp up period
#     tinds.append(np.arange(itsteady,itmax+1)) # ramp up period

#     return tinds

def prep_for_calc(dosim, iy, ix, it):
    '''
    Get drifter positions and unwrap them in the x direction for the periodic domain.

    (unwrap((d.variables['xg'][10,:]*(2*pi)/256.)))*(256./(2*pi))
    '''

    # Figure out which simulations to use based on time indices


    d = netCDF.MFDataset('tracks/balanced/' + dosim + '/')




def run(docalc, dosim):
    '''
    Run FSLE or dispersion calculation for sub regions of drifters in idealized shelf simulations.
    Inputs:
     docalc     'fsle' or 'dispersion'
     dosim      which simulation to use:
                'shelfstrat_M2_5.00e-07_N2_1.00e-04_f_5.00e-05'
                'shelfstrat_M2_3.54e-07_N2_1.00e-04_f_5.00e-05'
                'shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04'
                'shelfstrat_M2_7.07e-07_N2_1.00e-04_f_1.00e-04'
                'shelfstrat_M2_4.47e-07_N2_1.00e-04_f_1.00e-04'
                'shelfstrat_M2_3.16e-07_N2_1.00e-04_f_1.00e-04'
    '''

    # Each entry in strips list contains
    strips = select_pts()

    # # tinds[0] has ramp up time indices and tinds[1] has steady time indices
    # tinds = select_times(dosim)

    # run through time index sets
    for it, tind in xrange(2):

        # tinds[0] has ramp up time indices and tinds[1] has steady time indices
        tinds = select_times(dosim)

        # Run through selections out of drifter initialization locations in strips
        for i, strip in enumerate(strips):

            prep_for_calc(dosim, tind, strip[i][0], strip[i][1])
        
        # Run calculation

        # Average

        # Save






# if __name__ == "__main__":
#     run()    
