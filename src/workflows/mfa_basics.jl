

"""
compile_julia_fwsim( net::Network , target_semus::Vector{SpeciesEMU} )
  returns FWSim_Julia

creates the Julia_FWSim object for a specific network and a set of species emus that should be simulated
"""
function compile_julia_fwsim( net::Network , target_semus::Vector{SpeciesEMU} )

    emudecomp_1 = MetabolicFluxAnalysis.compute_emu_decomp( net , target_semus )
    emu_graph_1 = MetabolicFluxAnalysis.create_emu_graph(emudecomp_1)

    fluxi_ast_1 = MetabolicFluxAnalysis.create_ast_from_emu_graph( emu_graph_1 )

    # get the target semus from the fluxi measurements..
    target_smrs = MetabolicFluxAnalysis.construct_all_mass_ratios_for_measurement( target_semus )
    #julia_fw_sim_1 = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1,fluxi_ast_1,target_smrs)
    julia_fw_sim_1 = MetabolicFluxAnalysis.ast_to_juliacode_01(net,fluxi_ast_1,target_smrs)

    return julia_fw_sim_1
end

"""
compile_julia_fwsim( net::Network , data::Vector{MassDistributionMeasurement} )
  returns FWSim_Julia

creates the Julia_FWSim object for a specific network and a set of measurements
"""
function compile_julia_fwsim( net::Network , data::Vector{MassDistributionMeasurement}  )
    target_semus = map( d -> d.sass , data)
    julia_fw_sim =  compile_julia_fwsim( net , target_semus )
end


"""
struct ILE_Evaluator

contains all the information to evaluate flux distributions with respect to
a given ILE

Functions to use with:
  eval_flux2mds( ilee::ILE_Evaluator , x::Array{Float64,2} )
  eval_flux2sse( ilee , x::Array{Float64,2} )
"""
struct ILE_Evaluator
    substrate_mds::Array{Float64,1}
    mdata::Array{Float64,2}
    mdata2fw::Vector{Int64} # f_fwsim(x)[i]  corresponds to mdata[i,:]
    f_fwsim::Any           # function (substr_mds::Array{Float64,1},fluxes::Array{Float64,2}) -> mds::Array{Float64,2}
    f_fwsim_armed::Any
    fwsim::FWSim
end

"""
eval_flux2mds( ilee::ILE_Evaluator , x::Array{Float64,2} )
  returns Array{Float64,2}

computes mds for given ILE_Evaluator and fluxes
"""
function eval_flux2mds( ilee::ILE_Evaluator , x::Array{Float64,2} )
    mds::Array{Float64,2} = NaN * ones(size(ilee.mdata,1),size(x,2))
    mds                   = ilee.f_fwsim(ilee.substrate_mds,x)
    #for zi=1:size(x,2)
    #    mds[:,zi] = ilee.f_fwsim(ilee.substrate_mds,x[:,zi])
    #end
    return mds
end

"""
eval_flux2mds( ilee::ILE_Evaluator , x::Array{Float64,2} )
  returns Array{Float64,2}

computes sses for given ILE_Evaluator and fluxes
"""
function eval_flux2sse( ilee , x::Array{Float64,2} )
    mds = eval_flux2mds( ilee , x)
    res::Array{Float64,2} = broadcast( / , broadcast( - , mds[ilee.mdata2fw,:] , ilee.mdata[:,1] ) , ilee.mdata[:,2] )
    sse::Array{Float64,2} = sum( res.^2 , 1 )
    return sse
end

function eval_flux2sse( ilee , x::Array{Float64,1} )
    mds = eval_flux2mds( ilee , x[:,:])
    res::Array{Float64,2} = broadcast( / , broadcast( - , mds[ilee.mdata2fw,:] , ilee.mdata[:,1] ) , ilee.mdata[:,2] )
    sse::Array{Float64,2} = sum( res.^2 , 1 )[1]
    return sse
end


"""
create_julia_ile_evaluator( net::Network , substr_config::Vector{Dict{MetabolicFluxAnalysis.ConfiguredSEMU,Float64}} , data::Vector{MassDistributionMeasurement} )
  returns ILE_Evaluator

creates the ILE_Evaluator object which can be used to run the FW-Sim and perform
the fitting of fluxes against md data.
"""
function create_julia_ile_evaluator( net::Network , substr_config::Vector{Dict{MetabolicFluxAnalysis.ConfiguredSEMU,Float64}} , data::Vector{MassDistributionMeasurement} )
    julia_fw_sim = compile_julia_fwsim( net::Network , data::Vector{MassDistributionMeasurement}  )

    mr_positions::Dict{ASTEMUMassRatio,Int64} = Dict{ASTEMUMassRatio,Int64}()
    for zi in 1:length(julia_fw_sim.target_mass_ratios)
        mr_positions[julia_fw_sim.target_mass_ratios[zi]] = zi
    end

    mdata     = Array{Float64,2}(0,3)
    fw2mdata  = Vector{Int64}()

    for mdm in data
        semu = mdm.sass
        for zi=0:length(semu)
            mrd = ASTEMUMassRatio(semu,zi)
            if(!haskey(mr_positions,mrd))
                print("\n!! Something went wrong, no target semu for  measurement $(mrd)\n")
            else
                mdata    = [ mdata     ; mdm.data[zi][1:3]' ]
                fw2mdata = [ fw2mdata  ; mr_positions[mrd]  ]
            end
        end
    end

     # and collect data for semus:
     mds_a =  MetabolicFluxAnalysis.construct_mds( substr_config , julia_fw_sim.semu_mass_ratios )

    f_fwsim_single = eval(  parse( julia_fw_sim.code ) )
    f_fwsim_batch  = (xsubstr_mds::Array{Float64,1},xfluxes::Array{Float64,2}) -> eval_over_columns( (x_::Array{Float64,1}) -> f_fwsim_single(xsubstr_mds,x_) , size(mdata,1) , xfluxes )
    f_fwsim_batch_armed  = (xfluxes::Array{Float64,2}) -> eval_over_columns( (x_::Array{Float64,1}) -> f_fwsim_single(mds_a,x_) , size(mdata,1) , xfluxes )
    return ILE_Evaluator(mds_a,mdata,fw2mdata,f_fwsim_batch, f_fwsim_batch_armed, julia_fw_sim)
end

function eval_over_columns( f::Any , rows_result::Int , x::Array{Float64,2} )
    x_out::Array{Float64,2} = NaN * ones(rows_result,size(x,2))
    for zi=1:size(x,2)
        x_out[:,zi] = f(x[:,zi]) #ilee.f_fwsim(ilee.substrate_mds,x[:,zi])
    end
    return x_out
end


function fva(net::Network,lc::Vector{LinConstr})
    return fva(net::Network,lc::Vector{LinConstr},collect(1:length(net.fluxes)))
end

function fva(net::Network,lc::Vector{LinConstr},idx::Vector{Int64})
    print_output = length(idx)>10
    if(print_output); print("fva: "); end
    d  = length(net.fluxes)
    fp = get_flux_polyhedron(net,lc)
    lb = NaN * ones(length(idx))
    ub = NaN * ones(length(idx))
    for zi=1:length(idx)
        fba_result = solve_lp_minmax(fp, (1.0*( (1:d).== idx[zi] ) ) )
        lb[zi] = fba_result[1]
        ub[zi] = fba_result[2]
        if(print_output); if(mod(zi,ceil(length(idx)/40.))==0) print("."); end; end
    end
    if(print_output); print(" ..done!\n") end
    return (lb,ub)
end




"""
compile_opencl_fwsim( net::Network , target_semus::Vector{SpeciesEMU} )
  returns FWSim_OpenCL

creates the OpenCL_FWSim object for a specific network and a set of species emus that should be simulated
"""
function compile_opencl_fwsim( net::Network , target_semus::Vector{SpeciesEMU} )

    emudecomp_1 = MetabolicFluxAnalysis.compute_emu_decomp( net , target_semus )
    emu_graph_1 = MetabolicFluxAnalysis.create_emu_graph(emudecomp_1)

    fluxi_ast_1 = MetabolicFluxAnalysis.create_ast_from_emu_graph( emu_graph_1 )

    # get the target semus from the fluxi measurements..
    target_smrs     = MetabolicFluxAnalysis.construct_all_mass_ratios_for_measurement( target_semus )
    #julia_fw_sim_1 = MetabolicFluxAnalysis.ast_to_juliacode_01(fluxi1,fluxi_ast_1,target_smrs)
    opencl_fw_sim_1 = MetabolicFluxAnalysis.ast_to_opencl_01(net,fluxi_ast_1,target_smrs)

    return opencl_fw_sim_1
end

"""
compile_opencl_fwsim( net::Network , data::Vector{MassDistributionMeasurement} )
  returns FWSim_OpenCL

creates the OpenCL FWSim object for a specific network and a set of measurements
"""
function compile_opencl_fwsim( net::Network , data::Vector{MassDistributionMeasurement}  )
    target_semus = map( d -> d.sass , data)
    ocl_fw_sim =  compile_opencl_fwsim( net , target_semus )
end

function build_opencl_fwsim!( ocl_fwsim::MetabolicFluxAnalysis.FWSim_OpenCL , ocl_context::OCLContext)

    ocl_fwsim.ocl_context = ocl_context
    #p_ocl = cl.Program(ocl_context.ocl_ctx, source= ocl_fwsim.code ) |> OpenCL.cl.build!
    print("Create cl.Program\n")
    p_ocl = cl.Program(ocl_context.ocl_ctx, source= ocl_fwsim.code )
    print("Build! cl.Program\n")
    OpenCL.cl.build!(p_ocl)

    k_ocl = cl.Kernel(p_ocl, "fwsim")

    ocl_fwsim.p = p_ocl
    ocl_fwsim.k = k_ocl
end

function run_ocl_fluxes( ocl_fwsim::FWSim_OpenCL , semu_data_x::Array{Float64,1} , flux_data_x::Array{Float64,2} )
    n_fluxes = size(flux_data_x,2)
    data_out  = Array{Float32}( length(ocl_fwsim.target_mass_ratios) * n_fluxes )

    #semu_data = convert( Array{Float32,1} ,[0.00 ; 1.00; 1.00; 0.00; 0.00; 1.0 ;0.0;0.0;1.0;0.0;0.0;0;0;0;0;0] )
    #flux_data = convert( Array{Float32,1} ,[100.; 110 ; 50 ; 20 ; 20 ; 80 ; 0;0;0;0;0;0;0;0;0;0] )
    semu_data = convert( Array{Float32,1} , semu_data_x[:] )
    flux_data = convert( Array{Float32,1} , flux_data_x[:] )

    #device, ctx, queue = OpenCL.cl.create_compute_context()
    ocl_context = ocl_fwsim.ocl_context.value
    ctx         = ocl_context.ocl_ctx
    queue       = ocl_context.ocl_queue

    semu_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=semu_data)
    flux_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=flux_data)
    out_buff  = cl.Buffer(Float32, ctx, :w, length(data_out))
    k = ocl_fwsim.k.value #cl.Kernel(p, "fwsim")

    req_sx = ocl_fwsim.required_matrix_size_x
    req_sy = ocl_fwsim.required_matrix_size_y

    println("length data_out = $(length(data_out))")
    println("nfluxes = $(size(flux_data_x,2))")
    println("run config: $(req_sx),$(req_sy*size(flux_data,2)) ; $(req_sx),$(req_sy)")

    #queue(k, [ size(flux_data,2) * req_sx ; req_sy ], [req_sx;req_sy], semu_buff, flux_buff, out_buff)
    #queue(k, [ req_sx ; n_fluxes * req_sx ], [req_sx;req_sx], semu_buff, flux_buff, c_buff)
    queue(k, [ n_fluxes * req_sx ], [req_sx], semu_buff, flux_buff, out_buff)

    #cl.api.clGetKernelWorkGroupInfo(cl.CL_KERNEL_WORK_GROUP_SIZE)
    r = cl.read(ocl_context.ocl_queue, out_buff)
end

"""

"""
function run_ocl_fluxes( ocl_fwsim::FWSim_OpenCL , semu_buff::cl.Buffer , flux_data_x::Array{Float64,2} )
    n_fluxes = size(flux_data_x,2)
    data_out  = Array{Float32}( length(ocl_fwsim.target_mass_ratios) * n_fluxes )

    #semu_data = convert( Array{Float32,1} , semu_data[:] )
    flux_data = convert( Array{Float32,1} , flux_data_x[:] )

    ocl_context = ocl_fwsim.ocl_context.value

    #semu_buff = cl.Buffer(Float32, ocl_context.ctx, (:r, :copy), hostbuf=semu_data)
    flux_buff = cl.Buffer(Float32, ocl_context.ocl_ctx, (:r, :copy), hostbuf=flux_data)
    c_buff = cl.Buffer(Float32, ocl_context.ocl_ctx, :w, length(data_out))
    k = ocl_fwsim.k.value #cl.Kernel(p, "fwsim")

    req_sx = ocl_fwsim.required_matrix_size_x
    req_sy = ocl_fwsim.required_matrix_size_y
    #queue(k, [ req_sx ; n_fluxes * req_sx ], [req_sx;req_sx], semu_buff, flux_buff, c_buff)
    queue(k, [ n_fluxes * req_sx ], [req_sx], semu_buff, flux_buff, c_buff)
    r = cl.read(ocl_context.ocl_queue, c_buff)
end



"""
create_opencl_ile_evaluator( net::Network , substr_config::Vector{Dict{MetabolicFluxAnalysis.ConfiguredSEMU,Float64}} , data::Vector{MassDistributionMeasurement} )
  returns ILE_Evaluator

creates the ILE_Evaluator object which can be used to run the FW-Sim and perform
the fitting of fluxes against md data.
"""
function create_opencl_ile_evaluator( net::Network , substr_config::Vector{Dict{MetabolicFluxAnalysis.ConfiguredSEMU,Float64}} , data::Vector{MassDistributionMeasurement} , ocl_context::OCLContext  )
    ocl_fw_sim = compile_opencl_fwsim( net::Network , data::Vector{MassDistributionMeasurement}  )

    if(false)
        print("\nOCL KERNEL:\n")
        print(ocl_fw_sim.code)
        print("\n")
        return ocl_fw_sim
    end

    mr_positions::Dict{ASTEMUMassRatio,Int64} = Dict{ASTEMUMassRatio,Int64}()
    for zi in 1:length(ocl_fw_sim.target_mass_ratios)
        mr_positions[ocl_fw_sim.target_mass_ratios[zi]] = zi
    end

    mdata     = Array{Float64,2}(0,3)
    fw2mdata  = Vector{Int64}()

    for mdm in data
        semu = mdm.sass
        for zi=0:length(semu)
            mrd = ASTEMUMassRatio(semu,zi)
            if(!haskey(mr_positions,mrd))
                print("\n!! Something went wrong, no target semu for  measurement $(mrd)\n")
            else
                mdata    = [ mdata     ; mdm.data[zi][1:3]' ]
                fw2mdata = [ fw2mdata  ; mr_positions[mrd]  ]
            end
        end
    end

     # and collect data for semus:
     mds_a =  MetabolicFluxAnalysis.construct_mds( substr_config , ocl_fw_sim.semu_mass_ratios )

    #f_fwsim = eval(  parse( opencl_fw_sim.code ) )
    #p_ocl = cl.Program(ctx, source= opencl_fw_sim.code ) |> cl.build!
    #k_ocl = cl.Kernel(p, "fwsim")
    build_opencl_fwsim!( ocl_fw_sim::MetabolicFluxAnalysis.FWSim_OpenCL , ocl_context::OCLContext)

    # create the semu data buffer:
    semu_data = convert( Array{Float32,1} , mds_a[:] )
    semu_buff = cl.Buffer(Float32, ocl_context.ocl_ctx, (:r, :copy), hostbuf=semu_data)

    f_fwsim       = (xsubstr_mds::Array{Float64,1},xfluxes::Array{Float64,2}) -> run_ocl_fluxes(ocl_fw_sim,xsubstr_mds,xfluxes)
    f_fwsim_armed = (xfluxes::Array{Float64,2}) -> run_ocl_fluxes(ocl_fw_sim,semu_buff,xfluxes)
    #f_fwsim_batch  = (xsubstr_mds::Array{Float64,1},xfluxes::Array{Float64,2}) -> eval_over_columns( (x_::Array{Float64,1}) -> f_fwsim_single(xsubstr_mds,x_) , size(mdata,1) , xfluxes )

    return ILE_Evaluator(mds_a,mdata,fw2mdata,f_fwsim,f_fwsim_armed , ocl_fw_sim)
end
