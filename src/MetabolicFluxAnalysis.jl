module MetabolicFluxAnalysis


# package code goes here
export LC , Flux , Species , Network , add, parseLC , parse_linconstr , isless , compute_emu_decomp, EMURxn , getN
export SAS, RAS , PAS , SpeciesEMU  ,  constructSpeciesEMU , SpeciesEMURxn , create_emu_graph
export get_positions_of_semu
export parseFluxi , create_ast_from_emu_graph , parseFlux
export ASTEMUSim ,  ASTEMUSimStep , ASTMatrix , ASTNumNode , ASTEMUMassRatio , ASTFluxRate , ASTSymmOp , ASTRealValue

export construct_mds , construct_all_mass_ratios_for_measurement , get_configured_reactands

export lc2polyhedron , get_flux_polyhedron
export remove_redundant_constraints , mve_run , mve_solve , solve_chebyshev

# basic workflows:
export compile_julia_fwsim , evaluate_flux2mds , evaluate_flux2sse , create_julia_ile_evaluator
export ast_to_opencl_01

# io
export export_json , parse_ftbl_data , ftbl_to_fluxi , import_cobra

# wx parametrization:
export WXTransform , is_wx_flux , get_wx_flux_pair ,  get_wxflux_net , get_wxflux_xch , get_optimal_wx_parametrization

include("polyhedra/polyhedra.jl")
include("base/lincomb.jl")
include("base/isotopomer.jl")
include("base/networks.jl")
include("fluxi/fluxiparser.jl")
include("fluxi/emudecomp.jl")
include("fluxi/emugraph.jl")
include("fluxi/fluxi_ast.jl")

abstract type FWSim
end


include("fluxi/fluxi_julia_backend.jl")
include("fluxi/fluxi_opencl_backend.jl")

include("polyhedra/mve_solver.jl")
include("polyhedra/samp_unif.jl")
include("polyhedra/polytope_volume.jl")

include("workflows/mfa_basics.jl")
include("workflows/parametrization_basics.jl")

include("polyhedra/flux_parametrizations_01.jl")


include("io/json_io.jl")
include("io/import_ftbl.jl")
include("io/import_cobra.jl")

end # module
