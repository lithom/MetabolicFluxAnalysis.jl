
using JSON

"""
export_json(fd::FluxiData)
  returns String

creates the json representation of the fluxi data
"""
function export_json( fd::FluxiData )
    rxns = map( x::Flux -> export_json(x) , map( xs -> fd.net.fluxes_db[xs] ,fd.net.fluxes ) )
    species = map( x::Species -> export_json(x) , map( xs -> fd.net.species_db[xs] , fd.net.species ) )

    substrate_configs = Dict{String,Any}()
    for sci in collect(keys(fd.substrate_configs))
        substrate_configs[sci] = export_json(fd.substrate_configs[sci])
    end

    measurements = Dict{String,Any}()
    for msi in collect(keys(fd.measurements))
        measurements[msi] = export_json(fd.measurements[msi])
    end

    linconstraints = Dict{String,Any}()
    for lci in keys(fd.linear_constraints)
        linconstraints[lci] = export_json(fd.linear_constraints[lci])
    end

    json_buffer::IOBuffer = IOBuffer()

    fd = Dict{String,Any}()
    fd["reactions"] = rxns
    fd["metabolites"] = species
    fd["substrate_configs"] = substrate_configs
    fd["measurements"] = measurements
    fd["linear_constraints"] = linconstraints

    JSON.print(json_buffer, fd )

    return String(take!(json_buffer))
end

function export_json( f::Flux )
    fd = Dict{String,Any}()

    fd["name"] = f.id
    fd["id"]   = f.id

    if(!is_sx_flux(f))
        fd["lower_bound"] = 0.0
        fd["upper_bound"] = 1000.0

        fd["metabolites"] = add( f.P , mult( f.E , -1. ) )
        fd["metabolites_E"] = f.E
        fd["metabolites_P"] = f.P

        #fd["transition"] = fd.Transition
        # add atom transition data:
        fd["atomtransition"] = write_atom_transition(f)
    else
        fd["lower_bound"] = -1000.0
        fd["upper_bound"] =  1000.0

        if(haskey(f.E,"*"))
            fd["metabolites"] = f.P
            fd["metabolites_E"] = Dict{String,Any}()
            fd["metabolites_P"] = f.P
        else
            if(!haskey(f.P,"*"))
                error("Strange flux encountered: $(f)")
            end
            fd["metabolites"] = f.E
            fd["metabolites_E"] = Dict{String,Any}()
            fd["metabolites_P"] = f.E
        end
    end
    return fd
end


function export_json( s::Species )
    js = Dict{String,Any}()

    js["name"] = s.id
    js["id"]   = s.id
    js["notes"] = Dict{String,Any}( "CARBONS" => s.na )

    return js
end

"""
import_json_to_fluxi(str_json::String)
  returns FluxiData

creates the fluxi data object from the fluxi json representation
"""
function import_json_to_fluxi( str_json::String )
    imp = import_json(str_json)
    linear_constraints_sets = imp[7]
    linear_constraints_sets["rxn_bounds"] = imp[2]
    return FluxiData(imp[1],imp[5],imp[6],imp[7])
end

function import_json( str_json::String )
    json_d = JSON.parse(str_json)
    import_json(json_d)
end

function import_json( json_d::Dict{String,Any} )
    import_s = import_json_species(json_d["metabolites"])
    import_f = import_json_fluxes(json_d["reactions"])

    import_measurements      = import_json_measurements(json_d["measurements"])
    import_substrate_configs = import_json_substrconfig(json_d["substrate_configs"])
    import_linconstraints    = import_json_linconstraints(json_d["linear_constraints"])

    return ( Network(import_f[1],import_s[1], import_f[2],import_s[2]), import_f[4],import_s,import_f,import_measurements,import_substrate_configs,import_linconstraints)
end

function import_json_species( species_list::Vector{Any})
    map_id2speciesdata = Dict{String,Dict}()
    map_id2species     = Dict{String,Species}()
    list_species       = Vector{String}()

    for dsi in species_list
        id        = ""
        num_atoms = -1

        if( !isa(dsi,Dict) )
            error("unexpected species.. no dict?")
        end
        id = dsi["id"]
        map_id2speciesdata[id] = dsi

        if(haskey(dsi,"notes"))
            if( haskey(dsi["notes"],"CARBONS") )
                if(isa(dsi["notes"]["CARBONS"][1],String))
                    num_atoms = parse(Float64,dsi["notes"]["CARBONS"][1])
                end
                if(isa(dsi["notes"]["CARBONS"][1],Real))
                    num_atoms = Int(dsi["notes"]["CARBONS"][1])
                end
            end
        end
        #println(num_atoms)
        s1 = Species(id,Integer(round(num_atoms)))
        push!(list_species,id)
        map_id2species[id] = s1
    end

    return (list_species,map_id2species,map_id2speciesdata)
end

function import_json_fluxes( fluxes_list::Vector{Any})
    map_id2fluxdata   = Dict{String,Dict}()
    map_id2fluxes     = Dict{String,Flux}()
    list_fluxes       = Vector{String}()

    lin_constraints   = Vector{LinConstr}()

    for dfi in fluxes_list
        stoich::Dict{String,Float64} = Dict{String,Float64}()
        for ki in keys(dfi["metabolites"])
            stoich[ki] = dfi["metabolites"][ki]
        end
        if( dfi["lower_bound"] < 0 )
            #split into two..
            id = dfi["id"]
            id_fw = string(id,"_fw")
            id_bw = string(id,"_bw")

            fi_fw = Flux(id, stoich )
            fi_bw = Flux(id, mult(stoich,-1.) )

            map_id2fluxdata[id_fw] = dfi
            map_id2fluxdata[id_bw] = dfi
            push!(list_fluxes,id_fw)
            push!(list_fluxes,id_bw)
            map_id2fluxes[id_fw] = fi_fw
            map_id2fluxes[id_bw] = fi_bw

            push!(lin_constraints,parse_linconstr("$(id_fw) < $(dfi["upper_bound"])") )
            push!(lin_constraints,parse_linconstr("$(id_bw) < $(-dfi["lower_bound"])") )
        else
            #one flux
            id = dfi["id"]
            map_id2fluxdata[id] = dfi
            fi = Flux(id, stoich )
            push!(list_fluxes,id)
            map_id2fluxes[id] = fi

            push!(lin_constraints,parse_linconstr("$(id) < $(dfi["upper_bound"])") )
            if(dfi["lower_bound"]>0)
                push!(lin_constraints,parse_linconstr("$(id) > $(dfi["lower_bound"])") )
            end
        end
    end
    return (list_fluxes,map_id2fluxes,map_id2fluxdata,lin_constraints)
end


function export_json( substr_config:: Vector{ Dict{MetabolicFluxAnalysis.ConfiguredSEMU,Float64} } )
    #fd = Dict{String,Any}()
    substrates = Vector{Dict{String,Any}}()
    for conf_s in substr_config
        d_s = Dict{String,Any}()
        d_s["metabolite"] = collect( collect(keys(conf_s))[1].semu)[1].species
        d_s["pos"]        = get_positions_of_semu( collect(keys(conf_s))[1].semu )
        list_mix          = Vector{Dict{String,Any}}()
        for csr in keys( conf_s )
            d_mi = Dict{String,Any}()
            d_mi["conf"]  = csr.config
            d_mi["ratio"] = conf_s[csr]
            push!( list_mix , d_mi )
        end
        d_s["mix"] = list_mix
        push!( substrates , d_s)
    end
    return substrates
end

function array_to_intarray( va::Vector{Any})
    vi = Vector{Integer}(length(va))
    for zi in 1:length(va)
        vi[zi] = Integer(va[zi])
    end
    return vi
end

#function import_json_substrconfig( substr_conf::Dict{String,Vector{Dict{String,Any}}} )
function import_json_substrconfig( d_substr_conf::Dict{String,Any} )
    substrate_configs = Dict{String,Vector{Dict{MetabolicFluxAnalysis.ConfiguredSEMU,Float64}}}()
    for key in collect(keys(d_substr_conf))
        substr_i = d_substr_conf[key]
        substrate_config_i = Vector{Dict{MetabolicFluxAnalysis.ConfiguredSEMU,Float64}}()
        for si in substr_i
            semu_i = constructSpeciesEMU( si["metabolite"] , array_to_intarray( si["pos"] ) )
            configured_substrate_i = Dict{ConfiguredSEMU,Float64}()
            for dsc in si["mix"]
                csemu_i = ConfiguredSEMU(semu_i, array_to_intarray( dsc["conf"] ) )
                configured_substrate_i[csemu_i] = dsc["ratio"]
            end
            push!(substrate_config_i,configured_substrate_i)
        end
        substrate_configs[key] = substrate_config_i
    end
    return substrate_configs
end


function export_json( measurements_set::Vector{MetabolicFluxAnalysis.MassDistributionMeasurement} )
    measurements = Vector{Dict{String,Any}}()
    for mi in measurements_set
        dm = Dict{String,Any}()
        dm["metabolite"] = collect( mi.sass)[1].species
        dm["pos"]        = get_positions_of_semu( mi.sass )

        dm_data          = Vector{Dict{String,Any}}()
        for msi in keys(mi.data)
            d_data_i = Dict{String,Any}()
            d_data_i["mass"] = msi
            d_data_i["data"] = mi.data[msi]
            push!(dm_data,d_data_i)
        end
        dm["data"] = dm_data
        push!(measurements,dm)
    end
    return measurements
end

function import_json_measurements(d_measurements::Dict{String,Any})
    measurements = Dict{String,Vector{MetabolicFluxAnalysis.MassDistributionMeasurement}}()
    for key in collect(keys(d_measurements))
        meas_i = d_measurements[key]
        measurements_i = Vector{MetabolicFluxAnalysis.MassDistributionMeasurement}()
        for msi in meas_i
            semu_i = constructSpeciesEMU( msi["metabolite"] , array_to_intarray( msi["pos"] ) )
            mass_data_i = Dict{Int64,Vector{Float64}}()
            for mass_di in msi["data"]
                mass_data_i[Int64(mass_di["mass"])] = mass_di["data"]
            end
            mdi_i = MassDistributionMeasurement(semu_i,mass_data_i)
            push!(measurements_i,mdi_i)
        end
        measurements[key] = measurements_i
    end
    return measurements
end


function export_json( lcs::Vector{MetabolicFluxAnalysis.LinConstr} )
    vlcs = Vector{Dict{String,Any}}()
    for lci in lcs
        d_lci = Dict{String,Any}()
        d_lci["lhs"] = lci.lhs
        d_lci["rhs"] = lci.rhs
        d_lci["op"]  = lci.constr
        push!(vlcs,d_lci)
    end
    return vlcs
end

function export_json(lci::LC)
    return lci
end

function import_json_lc(d_lc::Dict{String,Any})
    lci = LC{String}()
    for ki in keys(d_lc)
        lci[ki] = Float64(d_lc[ki])
    end
    return lci
end

function import_json_linconstraints(d_lcs::Dict{String,Any})
    constraints = Dict{String,Vector{MetabolicFluxAnalysis.LinConstr}}()
    for key in collect(keys(d_lcs))
        lcs_i = d_lcs[key]
        linconstr_i = Vector{MetabolicFluxAnalysis.LinConstr}()

        for d_lci in lcs_i
            dlhs = import_json_lc(d_lci["lhs"])
            drhs = import_json_lc(d_lci["rhs"])
            d_op = 'x'
            if(isa(d_lci["op"],String)); d_op = d_lci["op"][1]; end
            if(isa(d_lci["op"],Char)); d_op   = d_lci["op"]; end
            lc_new = LinConstr( d_op , dlhs , drhs )
            push!(linconstr_i,lc_new)
        end
        constraints[key] = linconstr_i
    end
    return constraints
end
