using MetabolicFluxAnalysis
using CoordinateTransformations
using OpenCL

    # current tests..
    fluxi_1_f = open("testdata/test_01.fluxi")
    fluxi_1_s = readstring(fluxi_1_f)
    fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )

    fluxi1.net.fluxes_db[fluxi1.net.fluxes[1]].transition
    emudecomp_1 = MetabolicFluxAnalysis.compute_emu_decomp( fluxi1.net , [ fluxi1.measurements["MassDistribution_01"][1].sass ] )
    emu_graph_1 = MetabolicFluxAnalysis.create_emu_graph(emudecomp_1)



    configured_sas = MetabolicFluxAnalysis.get_configured_reactands(fluxi1,"SC_010")
    #fluxi_ast_1 = MetabolicFluxAnalysis.create_ast_from_emu_graph( emu_graph_1 , configured_sas )
    fluxi_ast_1 = MetabolicFluxAnalysis.create_ast_from_emu_graph( emu_graph_1 )

    unrolled_ast = Vector{Any}()
    MetabolicFluxAnalysis.ast_traverse_df!( fluxi_ast_1 , unrolled_ast)

    # get the target semus from the fluxi measurements..
    target_smrs = MetabolicFluxAnalysis.construct_all_mass_ratios_for_measurement(fluxi1,"MassDistribution_01")
    julia_fw_sim = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1.net,fluxi_ast_1,target_smrs)
    #julia_fw_sim_fast = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1.net,fluxi_ast_1,target_smrs;fast_lu=true)




    ocl_fw_sim = MetabolicFluxAnalysis.ast_to_opencl_01(fluxi1.net,fluxi_ast_1,target_smrs)

print("\n\n\nDA KERNEL..\n\n")
print(ocl_fw_sim.code)

kernel_str = ocl_fw_sim.code

#semu_data = convert( Array{Float32,1} , repmat([0.00 ; 1.00; 1.00; 0.00; 0.00; 1.0 ],1))
semu_data = convert( Array{Float32,1} , repmat([ 0.0 ; 1.0 ; 1.0 ; 0.0 ; 0.0 ; 1.0 ; 0.0 ; 0.0 ; 1.0 ; 0.0 ; 0.0 ],1))
flux_data = convert( Array{Float32,1} , repmat([100.; 110 ; 50 ; 20 ; 20 ; 80 ],4) )

a         = Array{Float32,1}(100)


ocl_devices = OpenCL.cl.devices()
ocl_device  = ocl_devices[2]
ocl_context =  MetabolicFluxAnalysis.OCLContext(ocl_device)

#device, ctx, queue = OpenCL.cl.create_compute_context()
device = ocl_context.ocl_device
ctx    = ocl_context.ocl_ctx
queue  = ocl_context.ocl_queue


semu_buff = OpenCL.cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=semu_data)
flux_buff = OpenCL.cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=flux_data)
c_buff = OpenCL.cl.Buffer(Float32, ctx, :w, length(a))
p = OpenCL.cl.Program(ctx, source=kernel_str)
OpenCL.cl.build!(p)

#p = cl.Program(ctx, source=kernel_02) |> cl.build!
print(kernel_str)
#cl.CL_device_info(1)

k = cl.Kernel(p, "fwsim")
#queue(k, [8;8], [8;8], semu_buff, flux_buff, c_buff)
queue(k, [8,8], [8,8], semu_buff, flux_buff, c_buff)

queue(k, [32,8], [8,8], semu_buff, flux_buff, c_buff)

#queue(k, [8;32], [8;8], semu_buff, flux_buff, c_buff)

#cl.api.clGetKernelWorkGroupInfo(cl.CL_KERNEL_WORK_GROUP_SIZE)
r = cl.read(queue, c_buff)

#cl.work_group_info( k , cl.CL_KERNEL_WORK_GROUP_SIZE ,device )
#cl.work_group_info( k , cl.CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE  ,device )


julia_fw_sim

##

test_A = [ -150.000   0.000   0.000  50.000   0.000
  0.000  -150.000   0.000   0.000  50.000
 20.000   0.000  -20.000   0.000   0.000
 110.000  20.000   0.000  -130.000   0.000
  0.000  110.000  20.000   0.000  -130.000]
test_B = [  0.000  -100.000
 -100.000   0.000
  0.000   0.000
  0.000   0.000
  0.000   0.000]
test_Y = [  1.000   0.000
  0.000   1.000]

test_X = test_A \ (test_B*test_Y)

lu_result = lu(test_A)

 (test_B*test_Y)

after_l_solve = lu_result[1] \ (test_B*test_Y)

after_u_solve = lu_result[2] \ (after_l_solve)

##

mA = [ -150.000  50.000
 110.000  -130.000]

mB = [ -100.000   0.000
  0.000  -20.000]

mY = [  0.000   1.000   0.000
  0.000   0.000   0.000]

mX2 = mA \ (mB * mY)
