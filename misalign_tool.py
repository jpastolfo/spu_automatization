#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 14:52:02 2022

@author: joao.astolfo
"""
# Packages:

import time
import numpy
import Shadow
import sys                                                                                                                                                                                          
sys.path.insert(0, '/home/ABTLUS/joao.astolfo/miniconda3/envs/name-env/lib/python3.7/site-packages/oasys_srw/')   
from optlnls.shadow import calc_und_flux
from run_SAPUCAIA import run_SAPUCAIA
from optlnls.source import get_k, und_source
from optlnls.importing import read_shadow_beam
from optlnls.math import get_fwhm
from optlnls.plot import plot_beam
from funcs import run_first_beam, mirror
import glob
import matplotlib.transforms as trans
# =============================================================================
# Fixed Parameters:
# =============================================================================

eBeamEnergy = 3.0 # storage ring energy [GeV]

e_spread = 0.00085 # energy Spread

per = 22.0e-3 # undulator period [m]

L = 2.244 # undulator length [m]

k_max = 1.71 # maximum deflection parameter

accept_hor = 60e-06 # beamline acceptance. Also used for calculating the undulator flux [rad]

accept_ver = 60e-06 # beamline acceptance. Also used for calculating the undulator flux [rad]  

source_beam = [64.7e-6, 5.1e-6, 3.8e-6, 1.4e-6] # electron beam parameters [Sx [m], Sy [m], Sx' [rad], Sy'[rad]]

# =============================================================================
# User Defined parameters:
# =============================================================================

energy_array = numpy.linspace(6000, 18000, 13) # photon beam energy. Also used for calculating h, K, B [eV]

atomic_plane = 'Si111' # crystal type. 'Si111' or 'Si311'

n_rays = 100000 # number of rays for shadow

nbins = 151 # bins for the undulator flux

plot_figure = True # plots beam size and divergence

show_plots = True # plots flux graphs

save_txt = True # saves a .txt file with beamline flux

open_tool = True

run = False

# =============================================================================
# Automatic Calculated parameters:
# =============================================================================

e_beam_size_X, e_beam_size_Z, e_beam_div_X, e_beam_div_Z = source_beam  

betaX = e_beam_size_X/e_beam_div_X

betaZ = e_beam_size_Z/e_beam_div_Z

emmX = e_beam_size_X*e_beam_div_X

emmZ = e_beam_size_Z*e_beam_div_Z

# =============================================================================
# Run:
# =============================================================================

rx_points = 6
ry_points = 6
rz_points = 6

rx_range = numpy.linspace(0, 450, rx_points)
ry_range = numpy.linspace(0, 90, ry_points)
rz_range = numpy.linspace(0, 85000, rz_points)

a = numpy.zeros((rx_points, ry_points, rz_points, 151, 151))

if run:

    t0 = time.time()
    
    flux_list = []; energy_list = []; resolution_list = [];
    
    delta_E = (11/13)*(15000/1000) - (55/13) + 3
    
    h, K, B = get_k(per, 'max', 15000, k_max)  
    
    size_X, div_X = und_source(emmX, betaX, e_spread, L, per, 15000, h)
    size_Z, div_Z = und_source(emmZ, betaZ, e_spread, L, per, 15000, h)
    
    beam = run_first_beam(n_rays=n_rays, energy=15000, delta_E=delta_E, size_x=size_X/1000, size_z=size_Z/1000, div_x=1e-06*div_X, div_z=1e-06*div_Z,
                        atomic_plane=atomic_plane, accept_hor=accept_hor, accept_ver=accept_ver)
    
    for i in range(0, len(rx_range)):
        for j in range(0, len(ry_range)):
            for k in range(0, len(rz_range)):
            
                copy_beam = beam.duplicate()
                
                new_beam = mirror(copy_beam, Tx_um=0, Ty_um=0, Tz_um=0, Rx_urad=rx_range[i], Ry_urad=ry_range[j], Rz_urad=rz_range[k])
                    
                beam2D_size = read_shadow_beam(new_beam, x_column_index=3, y_column_index=1, nbins_x=150, nbins_y=150, nolost=1, 
                                               ref=23, zeroPadding=0, gaussian_filter=0)
                
                
                a[i][j] = beam2D_size
                
                new = numpy.zeros((1, 151))
                new[0,0] = i
                new[0,1] = j
                new[0,2] = k
                new = numpy.append(new, beam2D_size, axis=0)
                
                #header = '{0:^5s}\n{1:^5s}'.format(str(i), str(j))
                numpy.savetxt('general-'+str(i)+'-'+str(j)+'-'+str(k)+'.txt', new)

    print('Calculation time = %.2f s' %(time.time()-t0))

if open_tool:        

    from matplotlib.widgets import Slider, Button
    from matplotlib import pyplot as plt
    
    filelist = glob.glob('general*.txt')
    
    for file in filelist:
        txt = numpy.genfromtxt(file)
        a[int(txt[0,0])][int(txt[0,1])][int(txt[0,2])] = txt[1:,:]
    
    # Define initial parameters
    init_rx = 0
    init_ry = 0
    init_rz = 0
    
    # Create the figure and the line that we will manipulate
    fig, ax = plt.subplots()
    mtx = numpy.transpose(numpy.transpose(a[0,0,0]))
    img = plt.pcolormesh(mtx[1:,0]*1e3, mtx[0,1:]*1e3, mtx[1:,1:])
    #width, height = fig.get_size_inches()       
    #fig.set_size_inches(width, width)
    
    # adjust the main plot to make room for the sliders
    plt.subplots_adjust(bottom=0.3)
    
    axtx = plt.axes([0.2, 0.2, 0.3, 0.03])
    tx_slider = Slider(
        ax=axtx,
        label=r'Tx [$\mu$m]',
        valmin=0,
        valmax=450,
        valstep=450/rx_points,
        valinit=init_rx,
    )
    
    axty = plt.axes([0.2, 0.15, 0.3, 0.03])
    ty_slider = Slider(
        ax=axty,
        label=r'Ty [$\mu$m]',
        valmin=0,
        valmax=90,
        valinit=init_ry,
        valstep=90/ry_points
    )
    
    axtz = plt.axes([0.2, 0.1, 0.3, 0.03])
    tz_slider = Slider(
        ax=axtz,
        label=r'Tz [$\mu$m]',
        valmin=0,
        valmax=90,
        valinit=init_ry,
        valstep=90/ry_points
    )
    
    axrx = plt.axes([0.65, 0.2, 0.3, 0.03])
    rx_slider = Slider(
        ax=axrx,
        label=r'Rx [$\mu$rad]',
        valmin=0,
        valmax=450,
        valstep=450/rx_points,
        valinit=init_rx,
    )
    
    axry = plt.axes([0.65, 0.15, 0.3, 0.03])
    ry_slider = Slider(
        ax=axry,
        label=r'Ry [$\mu$rad]',
        valmin=0,
        valmax=90,
        valinit=init_ry,
        valstep=90/ry_points
    )
    
    axrz = plt.axes([0.65, 0.1, 0.3, 0.03])
    rz_slider = Slider(
        ax=axrz,
        label=r'Rz [$\mu$rad]',
        valmin=0,
        valmax=85000,
        valinit=init_rz,
        valstep=85000/rz_points
    )
    
    # The function to be called anytime a slider's value changes
    def update(val):
        mtx = numpy.transpose(a[int(rx_points*rx_slider.val/450)][int(ry_points*ry_slider.val/90)][int(rz_points*rz_slider.val/85000)])
        ax.set_xbound(mtx[0,1:][0]*1e3, mtx[0,1:][-1]*1e3)
        ax.set_ybound(mtx[1:,0][0]*1e3, mtx[1:,0][-1]*1e3)
        img.set_transform(trans.Affine2D().translate(1e3*(mtx[0,1:][-1] + mtx[0,1:][0])/2, 1e3*(mtx[1:,0][-1] + mtx[1:,0][0])/2) + ax.transData)
        img.set_array(mtx[1:,1:])
        #img.set_transform(trans.Affine2D().translate(0, 900) + ax.transData)
        fig.canvas.draw_idle()
    
    # register the update function with each slider
    rx_slider.on_changed(update)
    ry_slider.on_changed(update)
    rz_slider.on_changed(update)
    
    # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', hovercolor='0.975')
    
    def reset(event):
        rx_slider.reset()
        ry_slider.reset()
        rz_slider.reset()
    button.on_clicked(reset)
    
    plt.show()