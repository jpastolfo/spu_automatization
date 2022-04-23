#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:56:30 2022

@author: humberto.junior
"""

# Packages:
    
from run_SPU import run_SPU
from optlnls.importing import read_srw_int, read_srw_wfr
from optlnls.plot import plot_beam


if(0):
    
    # Set Parameters:
        
    n_electrons = 100000 
    nx = 600 
    ny = 600 
    drift_after_M1 = 17.0 
    M1_error_filename = 'SPU_total_deformation_300mm_sh.dat' 
    multi_e = 1 
    plot_single_e = True
    
    
    # Run beamline on SRW (multi-e):
        
    run_SPU(n_electrons, nx, ny, drift_after_M1, M1_error_filename, multi_e, plot_single_e)


if(1):
    
    # Set Parameters:
        
    n_electrons = 100000 # Not used
    nx = 600 
    ny = 600 
    drift_after_M1 = 25.0 
    M1_error_filename = 'SPU_total_deformation_300mm_sh.dat' 
    multi_e = 0
    plot_single_e = True
    
    
    # Run beamline on SRW (single-e):
        
    wfr = run_SPU(n_electrons, nx, ny, drift_after_M1, M1_error_filename, multi_e, plot_single_e)    
    
    
    # Plot convolution:
    
    mtx_convol = read_srw_wfr(wfr, pol_to_extract=6, int_to_extract=1, unwrap=0) # int_to_extract=1 ---> means convolution
    
    prefix = M1_error_filename[:-6]
    
    pic_filename = prefix + 'convol.png'
    
    cmap= 'viridis' # 'jet'# 'plasma' # 'viridis'
    
    ur = 200e-03; rw = 200e-03;
    
    plot_beam(beam2D=mtx_convol, outfilename=pic_filename, cut=1, textA=1, textB=5, textC=7, textD=13, x_range=ur, y_range=ur,
    			  cmap=cmap, x_range_min=-rw, x_range_max=rw, y_range_min=-rw, y_range_max=rw,
    			  fitType=0, plot_title='Sapucaia: E = 8keV       Height Error')
                             

    
