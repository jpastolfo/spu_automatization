# Packages:

import os
import time
import copy

import sys
sys.path.insert(0, '/home/ABTLUS/humberto.junior/miniconda3/envs/py37/lib/python3.7/site-packages/oasys_srw/')  

try:
    from oasys_srw.srwlib import *
    from oasys_srw.uti_plot import *
except:
    from srwlib import *
    from uti_plot import *

import numpy


# Functions:
    
def from_shadow_to_matrix(file_path, unit_factor):
    
    """
    Imports a mirror height error file in shadow/dabam format and returns a 2D matrix (in meters) in which\
    the first row is the longitudinal coordinates, the first column is sagittal coordinates,\
    elements [1:,1:] are the corresponding height errors and element [0,0] is unused.\n
    :file_path: path to file including file name.
    :unit_factor: unit factor from the file units to meter, that is, 1e-3 if file is in [mm] for instance.
    
    """
    
    import numpy as np
    data = [] # array to store all elements in 1D list    
    with open(file_path, 'r') as datafile: # reads and append every element sequentially until readline() fails.
        while True:
            read_line = datafile.readline()
            if not read_line: break    
            for element in read_line.split():
                data.append(element)
    datafile.close()
    
    nl, nc = int(data[0]), int(data[1])
    matrix_data = np.zeros((nl+1, nc+1)) # allocate matrix
    matrix_data[0,:][1:] = np.array(data[2:nc+2], dtype='float')*unit_factor # [m] Longitudinal coordinates
    for line in range(nl): # associate elements in data array to the matrix lines 
        matrix_data[line+1,:] = np.array(data[1+(nc+1)*(line+1):1+(nc+1)*(line+2)], dtype='float')*unit_factor # [m]
    return matrix_data


def SRW_figure_error(file_name, unit_factor, angle_in, angle_out, orientation_x_or_y, crop=False, height_offset=False, L=5e3, W=5e3):
    
    """
    Returns an instance of srwlib.SRWLOptT() which simulates a mirror height error in SRW module.\
    To run this function, it is necessary to have srwlib able to import.\n
    :file_name: filename of file in shadow/dabam format
    :unit_factor: unit factor from the file units to meter, that is, 1e-3 if file is in [mm] for instance.
    :angle_in: incidence angle [rad] relative to mirror surface
    :angle_out: reflection angle [rad] relative to mirror surface
    :orientation_x_or_y: 'x' for horizontal or 'y' for vertical deflection
    :crop: (optional) if True, crops the matrix to new length L and width W and optionally offset data to make median value equal zero
    :height_offset: if True, a height offset will be added in all points so that the median value of central line is zero.
    :L: total mirror length in which the matrix must be cropped 
    :W: total mirror width in which the matrix must be cropped
    """
    
    height2D = from_shadow_to_matrix(file_name, unit_factor)
    if(crop):
        height2D_cropped = crop_height_error_matrix(height2D, L, W, height_offset)
        print('Actual L x W: {0:.3f} m x {1:.3f} m'.format(height2D_cropped[0,-1]-height2D_cropped[0,1], height2D_cropped[-1,0]-height2D_cropped[1,0]))
        return srwl_opt_setup_surf_height_2d(height2D_cropped, orientation_x_or_y, angle_in, angle_out)
    else:
        return srwl_opt_setup_surf_height_2d(height2D, orientation_x_or_y, angle_in, angle_out)
    

def run_SPU(n_electrons=1, nx=600, ny=600, drift_after_M1=17.0, M1_error_filename='', multi_e=0, plot_single_e=False):
    
    """
    Runs SAPUCAIA beamline on SRW. 
    
    """
    
    t0 = time.time()

    #if not srwl_uti_proc_is_master(): exit()
    
    ####################################################
    # LIGHT SOURCE
    
    part_beam = SRWLPartBeam()
    part_beam.Iavg               = 0.1
    part_beam.partStatMom1.x     = 0.0
    part_beam.partStatMom1.y     = 0.0
    part_beam.partStatMom1.z     = -1.21
    part_beam.partStatMom1.xp    = 0.0
    part_beam.partStatMom1.yp    = 0.0
    part_beam.partStatMom1.gamma = 5870.85355072162
    part_beam.arStatMom2[0]      = 4.18609e-09
    part_beam.arStatMom2[1]      = 0.0
    part_beam.arStatMom2[2]      = 1.444e-11
    part_beam.arStatMom2[3]      = 2.6010000000000002e-11
    part_beam.arStatMom2[4]      = 0.0
    part_beam.arStatMom2[5]      = 1.9599999999999997e-12
    part_beam.arStatMom2[10]     = 7.225e-07
    
    magnetic_fields = []
    magnetic_fields.append(SRWLMagFldH(1, 'v', 
                                        _B=0.7942323334589987, 
                                        _ph=0.0, 
                                        _s=-1, 
                                        _a=1.0))
    magnetic_structure = SRWLMagFldU(_arHarm=magnetic_fields, _per=0.022, _nPer=102.0)
    magnetic_field_container = SRWLMagFldC(_arMagFld=[magnetic_structure], 
                                            _arXc=array('d', [0.0]), 
                                            _arYc=array('d', [0.0]), 
                                            _arZc=array('d', [0.0]))
    
    mesh = SRWLRadMesh(_eStart=14999.997174425363,
                        _eFin  =14999.997174425363,
                        _ne    =1,
                        _xStart=-0.001,
                        _xFin  =0.001,
                        _nx    = nx, #600,
                        _yStart=-0.001,
                        _yFin  =0.001,
                        _ny    = ny, #600,
                        _zStart=26.0)
    
    stk = SRWLStokes()
    stk.allocate(1,600,600)
    stk.mesh = mesh
    
    wfr = SRWLWfr()
    wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
    wfr.mesh = mesh
    wfr.partBeam = part_beam
    wfr.unitElFld = 1
    
    initial_mesh = deepcopy(wfr.mesh)
    
    ####################################################
    # BEAMLINE
    
    srw_oe_array = []
    srw_pp_array = []
    
    oe_0=SRWLOptA(_shape='r',
                    _ap_or_ob='a',
                    _Dx=0.00156,
                    _Dy=0.00156,
                    _x=0.0,
                    _y=0.0)
    
    pp_oe_0 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]
    
    srw_oe_array.append(oe_0)
    srw_pp_array.append(pp_oe_0)
    
    
    drift_after_oe_1 = SRWLOptD(5.0)
    pp_drift_after_oe_1 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]
    
    srw_oe_array.append(drift_after_oe_1)
    srw_pp_array.append(pp_drift_after_oe_1)
    
    oe_2 = SRWLOptMirTor(_size_tang=0.3,
                          _size_sag=0.005,
                          _rt=7908.17941116253,
                          _rs=0.0968748022136628,
                          _ap_shape='r',
                          _sim_meth=2,
                          _treat_in_out=1,
                          _nvx=-0.9999938750062707,
                          _nvy=0,
                          _nvz=-0.003499992849008717,
                          _tvx=0.003499992849008717,
                          _tvy=0,
                          _x=0.0,
                          _y=0.0)
    oe_2.set_dim_sim_meth(_size_tang=0.3,
                          _size_sag=0.005,
                          _ap_shape='r',
                          _sim_meth=2,
                          _treat_in_out=1)
    oe_2.set_orient(_nvx=-0.9999938750062707,
                      _nvy=0,
                      _nvz=-0.003499992849008717,
                      _tvx=0.003499992849008717,
                      _tvy=0,
                      _x=0.0,
                      _y=0.0)
    
    
    pp_oe_2 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]
    
    srw_oe_array.append(oe_2)
    srw_pp_array.append(pp_oe_2)
    
    
    # Surface Error:
        
    prefix_error = ''
    
    if(M1_error_filename != ''):    
        M1_error = SRW_figure_error(file_name=M1_error_filename, unit_factor=1e-03, angle_in=3.5e-03, angle_out=3.5e-03, orientation_x_or_y='x')
        pp_M1_error = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]
    
        srw_oe_array.append(M1_error)
        srw_pp_array.append(pp_M1_error)
        
        prefix_error, ext = os.path.splitext(M1_error_filename)
        
    
    drift_after_oe_2 = SRWLOptD(drift_after_M1) #SRWLOptD(17.0)
    pp_drift_after_oe_2 = [0,0,1.0,1,0,1.0,2.0,1.0,2.0,0,0.0,0.0]
    
    srw_oe_array.append(drift_after_oe_2)
    srw_pp_array.append(pp_drift_after_oe_2)
    
    
    ####################################################
    # PROPAGATION
    
    optBL = SRWLOptC(srw_oe_array, srw_pp_array)
    
    
    ####################################################
    ### RUN SINGLE ELECTRON
    
    if not(multi_e):
        
        #### PERFORM ELECTRIC FIELD CALCULATION
        
        srwl.CalcElecFieldSR(wfr, 0, magnetic_field_container, [1,0.01,0.0,0.0,50000,1,0.0])
        
        
        xMesh = [wfr.mesh.xStart, wfr.mesh.xFin, wfr.mesh.nx]
        yMesh = [wfr.mesh.yStart, wfr.mesh.yFin, wfr.mesh.ny]
        intUnits = ['m', 'm', 'ph/s/.1%bw/mm^2']
        
        #### EXTRACT INTENSITY FROM ELECTRIC FIELD
        
        # HORIZONTAL POL
        arIx = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
        srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
        
        arPx = array('d', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
        srwl.CalcIntFromElecField(arPx, wfr, 0, 4, 3, wfr.mesh.eStart, 0, 0)        
           
        wfrc1 = deepcopy(wfr) #copy.deepcopy(wfr)
        srwl.PropagElecField(wfrc1, optBL)
        
        #### EXTRACT INTENSITY FROM   PROPAGATED  ELECTRIC FIELD
        
        # HORIZONTAL POL
        arIx = array('f', [0]*wfrc1.mesh.nx*wfrc1.mesh.ny) #"flat" array to take 2D intensity data
        srwl.CalcIntFromElecField(arIx, wfrc1, 6, 0, 3, wfrc1.mesh.eStart, 0, 0)
        
        arPx = array('d', [0]*wfrc1.mesh.nx*wfrc1.mesh.ny) #"flat" array to take 2D intensity data
        srwl.CalcIntFromElecField(arPx, wfrc1, 0, 4, 3, wfrc1.mesh.eStart, 0, 0)
        
        if(M1_error_filename != ''):
            single_e_filename = 'SAPUCAIA_'+str(round(drift_after_M1, 1)+31)+'m_'+ prefix_error +'_sE.dat'
        else:
            single_e_filename = 'SAPUCAIA_'+str(round(drift_after_M1, 1)+31)+'m_Ideal_Beam_sE.dat'    
            
        srwl_uti_save_intens_ascii(arIx, wfrc1.mesh, single_e_filename, 0)
        
        xMesh1 = [wfrc1.mesh.xStart, wfrc1.mesh.xFin, wfrc1.mesh.nx]
        yMesh1 = [wfrc1.mesh.yStart, wfrc1.mesh.yFin, wfrc1.mesh.ny]
        intUnits = ['m', 'm', 'ph/s/.1%bw/mm^2']

        print('total time = {0:.1f} s'.format(time.time()-t0))

        
        if(plot_single_e):
        
            from optlnls.importing import read_srw_wfr
            from optlnls.plot import plot_beam
            
            beam = read_srw_wfr(wfrc1, 6, 1)
    
            ur = 0
            wr = 0.2
            
            text = 'SPU: 15 keV - Sample     Ideal mirror'
            
            plot_beam(beam, units=2, unitFactor=1e3, textA=1, textB=5, textC=7, textD=13, plot_title=text, xlabel='Horizontal', ylabel='Vertical',
                      x_range=1, y_range=1, x_range_min=-wr, x_range_max=wr, y_range_min=-wr, y_range_max=wr, zero_pad_y=0, cmap='viridis')
            
            
        return wfrc1

        
    
    ####################################################
    ### RUN MULTI ELECTRON
    
    ### USE COMMAND:
    ### mpirun -n 8 python file_path.py
    ### changing "8" for number of processes and the python file name. 
    ### you will need openMPI and mpi4py
    
    if(multi_e):
        
        nMacroElec = n_electrons #total number of macro-electrons -> Maximum was done with 500
        nMacroElecAvgPerProc = 5 #number of macro-electrons / wavefront to average on worker processes before sending data to master (for parallel calculation only)
        nMacroElecSavePer = 10 #intermediate data saving periodicity (in macro-electrons)
        srCalcMeth = 1 #SR calculation method
        srCalcPrec = 0.01 #SR calculation rel. accuracy
        srSampFact = 0
        
        if(M1_error_filename != ''):
            multi_e_filename = 'SAPUCAIA_'+str(round(drift_after_M1, 1)+31)+'m_'+ prefix_error +'_mE.dat'
        else:
            multi_e_filename = 'SAPUCAIA_'+str(round(drift_after_M1, 1)+31)+'m_Ideal_Beam_mE.dat'                          

        radStokesProp = srwl_wfr_emit_prop_multi_e(part_beam,
                                                   magnetic_field_container,
                                                   wfr.mesh,
                                                   srCalcMeth,
                                                   srCalcPrec,
                                                   nMacroElec,
                                                   nMacroElecAvgPerProc,
                                                   nMacroElecSavePer,
                                                   multi_e_filename,
                                                   srSampFact,
                                                   optBL,
                                                   _char=0)
        
        print('total time = {0:.1f} s'.format(time.time()-t0))

    
    
    
    ####################################################
    # MULTI ELECTRON PROPAGATION
    
    # radStokesProp = srwl_wfr_emit_prop_multi_e(part_beam,
    #                                             magnetic_field_container,
    #                                             initial_mesh,
    #                                             1,
    #                                             0.01,
    #                                             500000,
    #                                             5,
    #                                             20,
    #                                             'output_srw_script_me.dat',
    #                                             1.0,
    #                                             optBL,
    #                                             _char=0)
