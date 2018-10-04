using MetabolicFluxAnalysis
using OpenCL

# current tests..
fluxi_1_f = open("testdata/test_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)
fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )


# create an evaluator for the forward simulation
ilee_jul = MetabolicFluxAnalysis.create_julia_ile_evaluator( fluxi1.net , fluxi1.substrate_configs["SC_010"] ,fluxi1.measurements["MassDistribution_01"] )

# create opencl evaluator for the fw simulation:
ocl_devices = OpenCL.cl.devices()
ocl_device  = ocl_devices[1]
ocl_context =  MetabolicFluxAnalysis.OCLContext(ocl_device)
ilee_ocl = MetabolicFluxAnalysis.create_opencl_ile_evaluator( fluxi1.net , fluxi1.substrate_configs["SC_010"] ,fluxi1.measurements["MassDistribution_01"] , ocl_context )
print("\nOpenCL Device: $(ocl_device)\n")

# create flux polytope:

fp1 = MetabolicFluxAnalysis.get_flux_polyhedron( fluxi1.net , fluxi1.linear_constraints["ConstraintsExchange"] , max_val=200 )

# parametrize the flux polytope
rfd_parametrization = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( fp1 ; remove_red_constr=1)

@assert MetabolicFluxAnalysis.isinside( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt[:,:] )


#                    v1  v2     v3   v4   v5   v6      Ain  Eout  Fout
reference_flux = [  100.; 110 ; 50 ; 20 ; 20 ; 80   ;  100 ; 60 ; 80 ]

ilee_jul.f_fwsim( ilee_jul.substrate_mds , reference_flux[:,[1]])
ilee_jul.f_fwsim( ilee_jul.substrate_mds , reference_flux[:,[1;1;1;1]])
ilee_ocl.f_fwsim( ilee_jul.substrate_mds , reference_flux[:,[1]])
ilee_ocl.f_fwsim( ilee_jul.substrate_mds , reference_flux[:,[1;1;1;1]])


flux_test_a = reference_flux[:,repmat([1],100000)]

tic(); results = ilee_jul.f_fwsim( ilee_jul.substrate_mds , flux_test_a ); toc()
tic(); results = ilee_ocl.f_fwsim( ilee_jul.substrate_mds , flux_test_a ); toc()




md_xx  = MetabolicFluxAnalysis.eval_flux2mds( ilee_jul , repmat(reference_flux[:,:],1,100) )
sse_xx = MetabolicFluxAnalysis.eval_flux2sse( ilee_jul , repmat(reference_flux[:,:],1,100) )


print(ilee_jul.fwsim.code)
print(ilee_ocl.fwsim.code)
