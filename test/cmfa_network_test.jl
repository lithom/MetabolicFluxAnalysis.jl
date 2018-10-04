


# load file
fluxi_1_f = open("testdata/cmfa_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)
fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )



fluxi_target_semus = reduce( (x,y)->[x;y] , map( x->x.sass , fluxi1.measurements["Data_Exp_01"] ) )
emudecomp_1 = MetabolicFluxAnalysis.compute_emu_decomp( fluxi1.net , fluxi_target_semus )
emu_graph_1 = MetabolicFluxAnalysis.create_emu_graph(emudecomp_1)

fluxi_ast_1 = MetabolicFluxAnalysis.create_ast_from_emu_graph( emu_graph_1 )

target_smrs = MetabolicFluxAnalysis.construct_all_mass_ratios_for_measurement(fluxi1,"Data_Exp_01")

julia_fw_sim_safe = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1.net,fluxi_ast_1,target_smrs;fast_lu=false)
julia_fw_sim_fast = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1.net,fluxi_ast_1,target_smrs;fast_lu=true)


jfw_safe = eval( parse( julia_fw_sim_safe.code ) )
jfw_fast = eval( parse( julia_fw_sim_fast.code ) )

mds_safe =  MetabolicFluxAnalysis.construct_mds(fluxi1.substrate_configs["SC_Conf_C1"],julia_fw_sim_safe.semu_mass_ratios)
mds_fast =  MetabolicFluxAnalysis.construct_mds(fluxi1.substrate_configs["SC_Conf_C1"],julia_fw_sim_fast.semu_mass_ratios)

print("\nInput Config:\n")
for zi=1:length(mds_safe)
    print("$(julia_fw_sim_safe.semu_mass_ratios[zi]) --> $(mds_safe[zi]) : $(mds_fast[zi]) \n")
end


## Create flux polyhedron:
constraints_config_01 = ["Constraints_FixedRatios","Constraints_BiomassRatios","Constraints_Measurements_01","Constraints_IrreversibilityPublication"]
constraints_01 = reduce( (x,y)->[x;y] , map( x -> fluxi1.linear_constraints[x] , constraints_config_01 ) )
fp1 = MetabolicFluxAnalysis.get_flux_polyhedron( fluxi1.net , constraints_01 , max_val=2000 )

# parametrize the flux polytope
rfd_parametrization = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( fp1 ; remove_red_constr=1)

@assert MetabolicFluxAnalysis.isinside( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt[:,:] )
@test MetabolicFluxAnalysis.isinside( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt[:,:] )
print("Flux Polytop has dimension $(size(rfd_parametrization.p_opt.G,2))")

xx = MetabolicFluxAnalysis.sample_unif_hr( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt , 5000 , 1000 );
xx_flux = MetabolicFluxAnalysis.eval_opt2raw(rfd_parametrization,xx)

results_safe = zeros(length(julia_fw_sim_safe.target_mass_ratios),5000)
results_fast = zeros(length(julia_fw_sim_fast.target_mass_ratios),5000)

println("\nBenchmark")
tic()
for zi=1:size(xx_flux,2)
    results_safe[:,zi] = jfw_safe(mds_safe,xx_flux[:,zi])
end
toc()
tic()
for zi=1:size(xx_flux,2)
    results_fast[:,zi] = jfw_fast(mds_fast,xx_flux[:,zi])
end
toc()

[ results_safe[:,1] results_fast[:,1] ]

maximum( abs.( results_safe[:]-results_fast[:] ) )

##
MetabolicFluxAnalysis.parametrization_wiechert_exchange( fluxi1.net , constraints_01 , 2000. ; sx_param="log" )
