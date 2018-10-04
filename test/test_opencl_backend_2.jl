using MetabolicFluxAnalysis
using OpenCL

fluxi_1_f = open("testdata/cmfa_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)
fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )



fluxi_target_semus = reduce( (x,y)->[x;y] , map( x->x.sass , fluxi1.measurements["Data_Exp_01"] ) )
emudecomp_1 = MetabolicFluxAnalysis.compute_emu_decomp( fluxi1.net , fluxi_target_semus )
emu_graph_1 = MetabolicFluxAnalysis.create_emu_graph(emudecomp_1)

fluxi_ast_1 = MetabolicFluxAnalysis.create_ast_from_emu_graph( emu_graph_1 )

target_smrs = MetabolicFluxAnalysis.construct_all_mass_ratios_for_measurement(fluxi1,"Data_Exp_01")
julia_fw_sim = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1.net,fluxi_ast_1,target_smrs)

## OpenCL:

ocl_devices = OpenCL.cl.devices()
ocl_device  = ocl_devices[2]
ocl_context =  MetabolicFluxAnalysis.OCLContext(ocl_device)
print("OpenCL Device: $(ocl_device)")

#semu_buff = cl.Buffer(Float32, ocl_context.ocl_ctx, (:r, :copy), hostbuf=semu_data)
#flux_buff = cl.Buffer(Float32, ocl_context.ocl_ctx, (:r, :copy), hostbuf=flux_data)
#c_buff = cl.Buffer(Float32, ocl_context.ocl_ctx, :w, length(a))
ocl_fw_sim = MetabolicFluxAnalysis.ast_to_opencl_01(fluxi1.net,fluxi_ast_1,target_smrs)

p = OpenCL.cl.Program(ocl_context.ocl_ctx, source=ocl_fw_sim.code)
OpenCL.cl.build!( p , options="-cl-opt-disable",raise = true)
k = cl.Kernel(p, "fwsim")
##

constraints_config_01 = ["Constraints_FixedRatios","Constraints_BiomassRatios","Constraints_Measurements_01","Constraints_IrreversibilityPublication"]
constraints_01 = reduce( (x,y)->[x;y] , map( x -> fluxi1.linear_constraints[x] , constraints_config_01 ) )
fp1 = MetabolicFluxAnalysis.get_flux_polyhedron( fluxi1.net , constraints_01 , max_val=2000 )
# parametrize the flux polytope
rfd_parametrization = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( fp1 ; remove_red_constr=1)
@assert MetabolicFluxAnalysis.isinside( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt[:,:] )
print("Flux Polytop has dimension $(size(rfd_parametrization.p_opt.G,2))")
xx = MetabolicFluxAnalysis.sample_unif_hr( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt , 5000 , 1000 );
xx_flux = MetabolicFluxAnalysis.eval_opt2raw(rfd_parametrization,xx)

ilee_jul = MetabolicFluxAnalysis.create_julia_ile_evaluator( fluxi1.net , fluxi1.substrate_configs["SC_Conf_C1"] ,fluxi1.measurements["Data_Exp_01"] )

tic()
@profile sse_xx = MetabolicFluxAnalysis.eval_flux2sse( ilee_jul , xx_flux )
toc()

tic()
semu_data = convert( Array{Float32,1} , ilee_jul.substrate_mds[:] )
flux_data = convert( Array{Float32,1} , xx_flux[:] )
data_out  = convert( Array{Float32,1} , zeros( size(xx_flux,2)*length(ilee_jul.fwsim.target_mass_ratios) ) )

#device, ctx, queue = OpenCL.cl.create_compute_context()
semu_buff = cl.Buffer(Float32, ocl_context.ocl_ctx, (:r, :copy), hostbuf=semu_data)
flux_buff = cl.Buffer(Float32, ocl_context.ocl_ctx, (:r, :copy), hostbuf=flux_data)
out_buff  = cl.Buffer(Float32, ocl_context.ocl_ctx, :w, length(data_out))
#k = ocl_fwsim.k.value #cl.Kernel(p, "fwsim")

req_sx = 256
req_sy = 256

println("length data_out = $(length(data_out))")
println("nfluxes = $(size(xx_flux,2))")
println("run config: $(req_sx),$(req_sy*size(xx_flux,2)) ; $(req_sx),$(req_sy)")

#queue(k, [ size(flux_data,2) * req_sx ; req_sy ], [req_sx;req_sy], semu_buff, flux_buff, out_buff)
ocl_context.ocl_queue(k, [ req_sx ; req_sy * size(xx_flux,2) ], [req_sx;req_sy], semu_buff, flux_buff, out_buff)

r = cl.read(queue, out_buff)


#sse_xx = MetabolicFluxAnalysis.eval_flux2sse( ilee_1 , xx_flux )
toc()
