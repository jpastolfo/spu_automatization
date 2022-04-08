def run_SAPUCAIA(n_rays=1000000, energy=12000, delta_E=8, size_x=0.065306, size_z=0.004812, div_x=9.626e-06, div_z=8.882e-06, atomic_plane='Si111', 
                 Tx_um=0, Ty_um=0, Tz_um=0, Rx_urad=0, Ry_urad=0, Rz_urad=0, accept_hor=60e-06, accept_ver=60e-06):
    import Shadow
    import numpy
    
    iwrite = 0

    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()
    
    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = accept_hor/2 #3e-05
    oe0.HDIV2 = accept_hor/2 #3e-05
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = n_rays #1000000
    oe0.PH1 = energy - (delta_E/2) #11996.0
    oe0.PH2 = energy + (delta_E/2) #12004.0
    oe0.SIGDIX = div_x #9.626e-06
    oe0.SIGDIZ = div_z #8.955e-06
    oe0.SIGMAX = size_x #0.064807
    oe0.SIGMAZ = size_z #0.006315
    oe0.VDIV1 = accept_ver/2 #3e-05
    oe0.VDIV2 = accept_ver/2 #3e-05
    
    oe1.DUMMY = 0.1
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.N_SCREEN = 1
    oe1.T_IMAGE = 3000.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 26000.0
    oe1.RX_SLIT = numpy.array([accept_hor * oe1.T_SOURCE, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.RZ_SLIT = numpy.array([accept_ver * oe1.T_SOURCE, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    
    oe2.DUMMY = 0.1
    oe2.FILE_REFL = bytes('/home/ABTLUS/joao.astolfo/Oasys/'+atomic_plane+'.dat', 'utf-8')
    oe2.FWRITE = 1
    oe2.F_CENTRAL = 1
    oe2.F_CRYSTAL = 1
    oe2.PHOT_CENT = energy #12000.0
    oe2.R_LAMBDA = 5000.0
    oe2.T_IMAGE = 50.0
    oe2.T_INCIDENCE = 80.5152566594
    oe2.T_REFLECTION = 80.5152566594
    oe2.T_SOURCE = 0.0
    
    oe3.ALPHA = 180.0
    oe3.DUMMY = 0.1
    oe3.FILE_REFL = bytes('/home/ABTLUS/joao.astolfo/Oasys/'+atomic_plane+'.dat', 'utf-8')
    oe3.FWRITE = 1
    oe3.F_CENTRAL = 1
    oe3.F_CRYSTAL = 1
    oe3.PHOT_CENT = energy #12000.0
    oe3.R_LAMBDA = 5000.0
    oe3.T_IMAGE = 1950.0
    oe3.T_INCIDENCE = 80.5152566594
    oe3.T_REFLECTION = 80.5152566594
    oe3.T_SOURCE = 0.0
    
    Tx = Tx_um/1000; Ty = Ty_um/1000; Tz = Tz_um/1000; Rx = numpy.rad2deg(1e-06*Rx_urad); Ry = numpy.rad2deg(1e-06*Ry_urad); Rz = numpy.rad2deg(1e-06*Rz_urad)
    
    oe4.ALPHA = 90.0
    oe4.DUMMY = 0.1
    oe4.FHIT_C = 1
    oe4.FILE_REFL = b'/home/ABTLUS/joao.astolfo/Oasys/Rh.dat'
    oe4.FMIRR = 3
    oe4.FWRITE = 1
    oe4.F_DEFAULT = 0
    oe4.F_MOVE = 1
    oe4.F_REFLEC = 1
    oe4.OFFX = -1*Ty
    oe4.OFFY = Tz
    oe4.OFFZ = -1*Tx
    oe4.RLEN1 = 150.0
    oe4.RLEN2 = 150.0
    oe4.RWIDX1 = 2.5
    oe4.RWIDX2 = 2.5
    oe4.SIMAG = 25000.0
    oe4.SSOUR = 31000.0
    oe4.THETA = 89.7994647717
    oe4.T_IMAGE = 25000.0
    oe4.T_INCIDENCE = 89.7994647717
    oe4.T_REFLECTION = 89.7994647717
    oe4.T_SOURCE = 0.0
    oe4.X_ROT = -1*Ry
    oe4.Y_ROT = Rz
    oe4.Z_ROT = -1*Rx
    
    #Run SHADOW to create the source
    if iwrite:
        oe0.write("start.00")
    beam.genSource(oe0)
    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #run optical element 1
    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")
    beam.traceOE(oe1,1)
    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")
    
    #run optical element 2
    print("    Running optical element: %d"%(2))
    if iwrite:
        oe2.write("start.02")
    beam.traceOE(oe2,2)
    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")

    #run optical element 3
    print("    Running optical element: %d"%(3))
    if iwrite:
        oe3.write("start.03")
    beam.traceOE(oe3,3)
    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")
        
    #run optical element 4
    print("    Running optical element: %d"%(4))
    if iwrite:
        oe4.write("start.04")
    beam.traceOE(oe4,4)
    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")
        
    return beam