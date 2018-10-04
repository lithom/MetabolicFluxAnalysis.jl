

#include("MetabolicFluxAnalysis.jl")
using MetabolicFluxAnalysis
reload("MetabolicFluxAnalysis")

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
@test 1 == 1

println(pwd())
path_project     = splitdir( pwd() )
path_to_datasets = joinpath( path_project[1] , "datasets")
println(path_to_datasets)


if(false)

    @testset "MFA.Basics_01" begin include("basic_tests_01.jl") end
    @testset "MFA.Basics_02" begin include("basic_tests_02.jl") end
    @testset "Fluxi.Basics_01" begin include("fluxi_tests_01.jl") end
    @testset "MFA.Polyhedra_01" begin include("tests_polyhedra.jl") end
    @testset "MFA.Sampling_01" begin include("test_samp.jl") end
    @testset "MFA.ImportCobra" begin include("test_cobra_import.jl") end
    @testset "Fluxi.AntoniewiczWorkflow" begin include("antoniewicz_toy_network_test.jl") end
    @testset "Fluxi.WXTransform_01" begin include("test_wx_transform.jl") end
    @testset "Fluxi.CMFA_Net_01" begin include("cmfa_network_test.jl") end
    @testset "Fluxi.ToyNetwork_01" begin include("fluxi_tests_01.jl") end
end

#@testset "Fluxi.AntoniewiczWorkflow2" begin include("antoniewicz_toy_network_test_02.jl") end
#@testset "Fluxi.FWSimTests_01" begin include("fw_sim_tests_01.jl") end
#@testset "MFA.ImportJSON" begin include("test_json_io.jl") end
@testset "Fluxi.OCL_01" begin include("test_opencl_backend.jl") end

# PARAMETRIATION FAILS SOMETIMES.. (Parametrization of polytope..) should try to make this work reliably..
@testset "MFA.ImportFTBL" begin include("test_ftbl_import.jl") end

# @testset "Fluxi.OpenCL_01" begin include("kernel_compilation_01.jl") end # does not work, not needed anymore..



#include("basic_tests_01.jl")
#include("basic_tests_02.jl")
#include("fluxi_tests_01.jl")
#include("antoniewicz_toy_network_test.jl")
