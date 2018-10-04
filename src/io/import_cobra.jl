
using MetabolicFluxAnalysis
using Requests
using MAT


function import_cobra( filename::String , modelname::String )

    mf = MAT.matopen( filename )
    mv = MAT.matread( filename )

    model     = mv[modelname]
    modelKeys = keys(model)

    species_1   = model["metNames"]
    fluxes_pre  = model["rxns"]

    S    = model["S"]
    revs = model["rev"]


    # create species db
    species        = Vector{String}(length(species_1))
    species_db     = Dict{String,Species}()
    for zi=1:length(species_1)
        species_db[species_1[zi]] = Species(species_1[zi],0)
        species[zi]               = species_1[zi]
    end

    # create fluxes
    fluxes      = Vector{String}()
    flux_db     = Dict{String,Flux}()

    for zi=1:length(fluxes_pre)

        idx_neg = find( S[:,zi] .< 0 )
        idx_pos = find( S[:,zi] .> 0 )
        e::LC = LC{String}(); for i_n in idx_neg; e[species[i_n]] = -S[i_n,zi]; end
        p::LC = LC{String}(); for i_p in idx_pos; p[species[i_p]] =  S[i_p,zi]; end
        if(isempty(e)); e["*"]=1.; end
        if(isempty(p)); p["*"]=1.; end

        if(revs[zi]>0)
            push!(fluxes,string(fluxes_pre[zi],"_fw"))
            push!(fluxes,string(fluxes_pre[zi],"_bw"))
            e2 = copy(p)
            p2 = copy(e)
            flux_db[string(fluxes_pre[zi],"_fw")] = Flux(fluxes_pre[zi],e,p,Dict{RAS,PAS}())
            flux_db[string(fluxes_pre[zi],"_bw")] = Flux(fluxes_pre[zi],e2,p2,Dict{RAS,PAS}())
        else
            push!(fluxes,fluxes_pre[zi])
            flux_db[fluxes_pre[zi]] = Flux(fluxes_pre[zi],e,p,Dict{RAS,PAS}())
        end
    end

    return Network(fluxes,species,flux_db,species_db)
end
