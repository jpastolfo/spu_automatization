# Packages:
import time
import numpy as np
import sys                                                                                    
sys.path.insert(0, '/home/ABTLUS/joao.astolfo/miniconda3/envs/name-env/lib/python3.7/site-packages/oasys_srw/')   
from optlnls.shadow import calc_und_flux
from run_SAPUCAIA import run_SAPUCAIA
from optlnls.source import get_k, und_source
from optlnls.importing import read_shadow_beam
from optlnls.math import get_fwhm
from optlnls.plot import plot_beam

# Fixed Parameters:
eBeamEnergy = 3.0 # storage ring energy [GeV]
e_spread = 0.00085 # energy Spread
per = 22.0e-3 # undulator period [m]
L = 2.244 # undulator length [m]  
k_max = 1.71 # maximum deflection parameter
accept_hor = 60e-06 # beamline acceptance. Also used for calculating the undulator flux [rad]
accept_ver = 60e-06 # beamline acceptance. Also used for calculating the undulator flux [rad]  
source_beam = [64.7e-6, 5.1e-6, 3.8e-6, 1.4e-6] # electron beam parameters [Sx [m], Sy [m], Sx' [rad], Sy'[rad]]

# User Defined parameters:
energy_array = np.linspace(6000, 18000, 13) # photon beam energy. Also used for calculating h, K, B [eV]
atomic_plane = 'Si111' # crystal type. 'Si111' or 'Si311'
n_rays = 1000000 # number of rays for shadow
nbins = 151 # bins for the undulator flux
plot_figure = True # plots beam size and divergence
show_plots = True # plots flux graphs
save_txt = True # saves a .txt file with beamline flux

# Automatic Calculated parameters:
e_beam_size_X, e_beam_size_Z, e_beam_div_X, e_beam_div_Z = source_beam  
betaX = e_beam_size_X/e_beam_div_X  
betaZ = e_beam_size_Z/e_beam_div_Z
emmX = e_beam_size_X*e_beam_div_X 
emmZ = e_beam_size_Z*e_beam_div_Z

# Run:
t0 = time.time()

flux_list = []; energy_list = []; resolution_list = [];

for energy in energy_array:
    h, K, B = get_k(per, 'max', energy, k_max)
    size_X, div_X = und_source(emmX, betaX, e_spread, L, per, energy, h)
    size_Z, div_Z = und_source(emmZ, betaZ, e_spread, L, per, energy, h)
    delta_E = (11/13)*(energy/1000) - (55/13) + 3
    
    beam = run_SAPUCAIA(n_rays=n_rays, energy=energy, delta_E=delta_E, size_x=size_X/1000, size_z=size_Z/1000, div_x=1e-06*div_X, div_z=1e-06*div_Z,
                        atomic_plane=atomic_plane, Tx_um=0, Ty_um=0, Tz_um=0, Rx_urad=0, Ry_urad=0, Rz_urad=0, accept_hor=accept_hor, accept_ver=accept_ver)

    outputs = calc_und_flux(beam=beam, nbins=nbins, eBeamEnergy=eBeamEnergy, eSpread=e_spread, current=0.1,
                            und_per=per, und_length=L, B=B, min_harmonic=1, max_harmonic=(h+2),
                            source_beam=source_beam, show_plots=show_plots, accept_hor=accept_hor, accept_ver=accept_ver)
    
    flux_s, flux_b = outputs['total flux at source'], outputs['total flux propagated']
    power_s, power_b = outputs['total power at source'], outputs['total power propagated']
    E_b, T_E = outputs['energy array'] , outputs['transmission array']
    
    if(plot_figure):
        beam2D_size = read_shadow_beam(beam, x_column_index=3, y_column_index=1, nbins_x=200, nbins_y=150, nolost=1, 
                                       ref=23, zeroPadding=0, gaussian_filter=0)
        
        beam2D_div = read_shadow_beam(beam, x_column_index=6, y_column_index=4, nbins_x=150, nbins_y=150, nolost=1, 
                                      ref=23, zeroPadding=0, gaussian_filter=0)
        
        if((energy/1000) < 10):
            title = 'SAP - %.1f keV            F: %.2E ph/s/100mA' %(energy/1000, flux_b)
        else:
            title = 'SAP - %.1f keV           F: %.2E ph/s/100mA' %(energy/1000, flux_b)
        
        plot_beam(beam2D_size, outfilename='SAP_beam_size_'+atomic_plane+'_'+str(int(energy))+'eV.png', 
                  cut=0, textA=1, textB=5, textC=6, textD=10, fitType=3, xlabel='Hor.', ylabel='Vert.', plot_title=title, unitFactor=1e3, 
                  x_range=1, y_range=1, x_range_min=-200e-3, x_range_max=200e-3, y_range_min=-200e-3, y_range_max=200e-3, cmap='viridis', 
                  integral=flux_b, zero_pad_x=1, zero_pad_y=4, export_slices=0)
        
        plot_beam(beam2D_div, outfilename='SAP_beam_div_'+atomic_plane+'_'+str(int(energy))+'eV.png', cut=0, 
                  textA=1, textB=5, textC=6, textD=10, fitType=3, xlabel='Hor.', ylabel='Vert.', plot_title=title, unitFactor=1e6, 
                  fwhm_threshold=0.5, x_range=1, y_range=1, x_range_min=-50e-6, x_range_max=50e-6, y_range_min=-50e-6, y_range_max=50e-6, 
                  cmap='viridis', integral=flux_b, units='$\mu$rad', zero_pad_x=1, zero_pad_y=1, export_slices=0)
    
    resol = get_fwhm(E_b, T_E)[0]/energy
    flux_list.append(flux_b)
    energy_list.append(energy)
    resolution_list.append(resol)
        
# Writing Flux .txt file:
if(save_txt):
    filename = 'SAP_flux_'+atomic_plane+'.txt'
    data = np.array([energy_list, flux_list, resolution_list])
    np.savetxt(filename, data.transpose(), '%.6E', header='Energy[eV]\tFlux[ph/s/100mA]\tResolution')

print('Total time = %.2f s' %(time.time()-t0))