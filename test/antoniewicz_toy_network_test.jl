using MetabolicFluxAnalysis
using CoordinateTransformations

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


    # get the target semus from the fluxi measurements..
    target_smrs = MetabolicFluxAnalysis.construct_all_mass_ratios_for_measurement(fluxi1,"MassDistribution_01")
    julia_fw_sim_safe = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1.net,fluxi_ast_1,target_smrs)
    julia_fw_sim_fast = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1.net,fluxi_ast_1,target_smrs;fast_lu=true)

    ocl_fw_sim_fast = MetabolicFluxAnalysis.ast_to_opencl_01(fluxi1.net,fluxi_ast_1,target_smrs)


    julia_fw_sim_1 = julia_fw_sim_safe

    println("\n\nJULIA CODE:\n\n")
    print(julia_fw_sim_safe.code)

    jfw1 = eval( parse( julia_fw_sim_safe.code ) )
    jfw1(rand(11),rand(9))


    mds_a =  MetabolicFluxAnalysis.construct_mds(fluxi1.substrate_configs["SC_010"],julia_fw_sim_1.semu_mass_ratios)
    print("\nInput Config:\n")
    for zi=1:length(mds_a)
        print("$(julia_fw_sim_1.semu_mass_ratios[zi]) --> $(mds_a[zi])\n")
    end

    # check if everythings right:
    #                            v1  v2    v3   v4   v5   v6      Ain  Eout  Fout
    results_test = jfw1(mds_a,[ 100.; 110 ; 50 ; 20 ; 20 ; 80   ;  100 ; 60 ; 80 ])

    #m_ref_01 = fluxi1.measurements["MassDistribution_01"]
    m1_ref = MetabolicFluxAnalysis.get_mass_distribution_measurement_value(fluxi1,"MassDistribution_01",target_smrs[1])
    m2_ref = MetabolicFluxAnalysis.get_mass_distribution_measurement_value(fluxi1,"MassDistribution_01",target_smrs[2])
    m3_ref = MetabolicFluxAnalysis.get_mass_distribution_measurement_value(fluxi1,"MassDistribution_01",target_smrs[3])
    m4_ref = MetabolicFluxAnalysis.get_mass_distribution_measurement_value(fluxi1,"MassDistribution_01",target_smrs[4])

    @test abs( results_test[1] - m1_ref[1,1] ) < 0.0001
    @test abs( results_test[2] - m2_ref[1,1] ) < 0.0001
    @test abs( results_test[3] - m3_ref[1,1] ) < 0.0001
    @test abs( results_test[4] - m4_ref[1,1] ) < 0.0001


    #reload("MetabolicFluxAnalysis")
    # now get and then parametrize the flux polytope
    fp = MetabolicFluxAnalysis.get_flux_polyhedron( fluxi1.net , fluxi1.linear_constraints["ConstraintsExchange"] , max_val=200 )

    fp_min = MetabolicFluxAnalysis.find_minimal_affine_subspace(fp)



    if(false)
        # benchmark: # currently this does not work nicely.. we need steady state flux distr. to get ok results..
        n_benchmark_runs = 200_000
        fluxes_i  = rand(9,n_benchmark_runs) + 0.02;
        #smrs_i    = rand(11,n_benchmark_runs) + 0.005;
        results_i = NaN * ones(4,n_benchmark_runs);
        tic()
            for zi=1:size(fluxes_i,2)
                results_i[:,zi] = jfw1(mds_a,fluxes_i[:,zi])
            end
        time_1 = toc()
        println("Eval Speed: $(n_benchmark_runs/time_1) / sec")
    end

# run this benchmark..
MetabolicFluxAnalysis.benchmark_matrix_solving_01()
# and run this..
MetabolicFluxAnalysis.test_lu_solving_01()
##


#wccs = weakly_connected_components( gg1.g )
