
##
reload("MetabolicFluxAnalysis")
using MetabolicFluxAnalysis


# load file
fluxi_1_f = open("datasets/complete_mfa/cmfa_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)
fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )

# create an evaluator for the forward simulation
ilee_1 = MetabolicFluxAnalysis.create_julia_ile_evaluator( fluxi1.net , fluxi1.substrate_configs["SC_Conf_C1"] ,fluxi1.measurements["Data_Exp_01"] )

ocl_devices = OpenCL.cl.devices()
ocl_device  = ocl_devices[3]
print("\nOpenCL Device: $(ocl_device)\n")
ilee_2 = MetabolicFluxAnalysis.create_opencl_ile_evaluator( fluxi1.net , fluxi1.substrate_configs["SC_Conf_C1"] ,fluxi1.measurements["Data_Exp_01"] , MetabolicFluxAnalysis.OCLContext(ocl_device))


f = open("cl_kernel_01.ocl","w")
write(f,ilee_2.code)

constraints_config_01 = ["Constraints_FixedRatios","Constraints_BiomassRatios","Constraints_Measurements_01","Constraints_IrreversibilityPublication"]
constraints_01 = reduce( (x,y)->[x;y] , map( x -> fluxi1.linear_constraints[x] , constraints_config_01 ) )
fp1 = MetabolicFluxAnalysis.get_flux_polyhedron( fluxi1.net , constraints_01 , max_val=2000 )

# parametrize the flux polytope
rfd_parametrization = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( fp1 ; remove_red_constr=1)

@assert MetabolicFluxAnalysis.isinside( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt[:,:] )
print("Flux Polytop has dimension $(size(rfd_parametrization.p_opt.G,2))")

xx = MetabolicFluxAnalysis.sample_unif_hr( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt , 5000 , 1000 );
xx_flux = MetabolicFluxAnalysis.eval_opt2raw(rfd_parametrization,xx)

tic()
@profile sse_xx = MetabolicFluxAnalysis.eval_flux2sse( ilee_1 , xx_flux )
toc()
