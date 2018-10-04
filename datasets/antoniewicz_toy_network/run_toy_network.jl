
##
reload("MetabolicFluxAnalysis")
using MetabolicFluxAnalysis




# load file
fluxi_1_f = open("datasets/antoniewicz_toy_network/toy_network_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)
fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )

# create an evaluator for the forward simulation
ilee_1 = MetabolicFluxAnalysis.create_julia_ile_evaluator( fluxi1.net , fluxi1.substrate_configs["SC_010"] ,fluxi1.measurements["MassDistribution_01"] )

# create flux polytope:

fp1 = MetabolicFluxAnalysis.get_flux_polyhedron( fluxi1.net , fluxi1.linear_constraints["ConstraintsExchange"] , max_val=200 )

# parametrize the flux polytope
rfd_parametrization = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( fp1 ; remove_red_constr=1)

@assert MetabolicFluxAnalysis.isinside( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt[:,:] )


#                    v1  v2     v3   v4   v5   v6      Ain  Eout  Fout
reference_flux = [  100.; 110 ; 50 ; 20 ; 20 ; 80   ;  100 ; 60 ; 80 ]

ilee_1.f_fwsim( ilee_1.substrate_mds , reference_flux[:,[1;1;1;1]])

md_xx  = MetabolicFluxAnalysis.eval_flux2mds( ilee_1 , repmat(reference_flux[:,:],1,100) )
sse_xx = MetabolicFluxAnalysis.eval_flux2sse( ilee_1 , repmat(reference_flux[:,:],1,100) )

# generate some unifrom samples over polytope:
xx = MetabolicFluxAnalysis.sample_unif_hr( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt , 5000 , 1000 );
xx_flux = MetabolicFluxAnalysis.eval_opt2raw(rfd_parametrization,xx)


sse_xx = MetabolicFluxAnalysis.eval_flux2sse( ilee_1 , xx_flux )
#sse_xx = MetabolicFluxAnalysis.eval_flux2sse( ilee_1 , reference_flux[:,:] )

using Plots
gr()
Plots.scatter( xx_flux[2,:],xx_flux[5,:] , zcolor=sse_xx[:],c=:lightrainbow,ms=4.0,marker=(stroke(0,:gray)),clims=(0,20))

#Plots.plot!()


using Plots
gr()
plot(xx_flux[2,:],xx_flux[5,:],sse_xx,nlevels=10)
ipX = ipY = linspace(-3, 3, 7)
ipZ = GR.interp2(xx_flux[2,:],xx_flux[5,:],sse_xx,ipX,ipY )





default(size=(500,300), leg=true)
y = rand(100)
scatter(y, z=abs(y-.5), m=10, c=:heat, lab="grad")
plot!(0:10:100,rand(11,4),lab="lines", palette=:grays)



# benchmark of unif. sampling:
tic(); xx = MetabolicFluxAnalysis.sample_unif_hr( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt , 1000 ,1000 ); toc();




using MetabolicFluxAnalysisWorkflows
f_cost = (x) -> MetabolicFluxAnalysis.eval_flux2sse( ilee_1 , MetabolicFluxAnalysis.eval_opt2raw(rfd_parametrization,x[:,:]) )[1]
opti_result = MetabolicFluxAnalysisWorkflows.find_optimum( rfd_parametrization.p_opt  , f_cost  , xx[:,1] )
opti_flux = MetabolicFluxAnalysis.eval_opt2raw( rfd_parametrization , opti_result[2] )
scatter!(opti_flux[2:2],opti_flux[5:5])
opti_result = MetabolicFluxAnalysisWorkflows.find_optimum( rfd_parametrization.p_opt  , f_cost  , xx )


#julia_fwsim_1 =  MetabolicFluxAnalysis.compile_julia_fwsim( fluxi1.net , fluxi1.measurements["MassDistribution_01"] )
#emudecomp_1 = MetabolicFluxAnalysis.compute_emu_decomp( fluxi1.net , [ fluxi1.measurements["MassDistribution_01"][1].sass ] )
#emu_graph_1 = MetabolicFluxAnalysis.create_emu_graph(emudecomp_1)
#fluxi_ast_1 = MetabolicFluxAnalysis.create_ast_from_emu_graph( emu_graph_1 )
#julia_fw_sim_1 = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1,fluxi_ast_1,target_smrs)


#str_json_01 = MetabolicFluxAnalysis.export_json(fluxi1)
