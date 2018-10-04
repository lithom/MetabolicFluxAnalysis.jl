import Base.isequal
import Base.hash
import Base.show
import Base.isless

"""
EMURxn has fields:
    flux::String, product_emu::ProductEMU, reactand_emus::Vector{ReactandEMU}

and represents an EMU Reaction
"""
struct EMURxn
    flux::String
    product_emu::ProductEMU
    reactand_emus::Vector{ReactandEMU}
    #transition::Dict{RAS,PAS}
end

"""
SpeciesEMURxn has fields:
    flux::String, product_emu::SpeciesEMU, reactand_emus::Vector{SpeciesEMU}

and represents an EMU Reaction on the species / SAS level
"""
struct SpeciesEMURxn
    flux::String
    product_emu::SpeciesEMU
    reactand_emus::Vector{SpeciesEMU}
    #transition::Dict{RAS,PAS}
end

# order: 1. product emu, 2. flux , this should be ok
function Base.isless(a::SpeciesEMURxn,b::SpeciesEMURxn)
    if(a.product_emu<b.product_emu) ; return true; end
    if(a.product_emu==b.product_emu)
        return a.flux < b.flux
    else
        return false
    end
end
function Base.isequal(a::SpeciesEMURxn,b::SpeciesEMURxn)
    if(a.flux < b.flux )
        return true
    end
    if(f.flux == b.flux)
        if( isequal(a.product_emu , b.product_emu) && a.flux ==b.flux )
            if(length(a.reactand_emus)!=lenght(b.reactand_emus))
                return false
            else
                rea = sort(a.reactand_emus)
                reb = sort(b.reactand_emus)
                equal = all( isequal.(rea,reb) )
                return equal
            end
        else
            return false
        end
    else
        return false
    end
    return false
end
function Base.hash(x::SpeciesEMURxn)
    hash( string("semurxn$(x.flux):$(hash(x.Set(product_emu))):$(hash(Set(x.educt_emus)))" )  )
end




function SpeciesEMURxn(r::EMURxn)
    SpeciesEMURxn( r.flux , pemu2semu( r.product_emu ) , map( x1 -> remu2semu(x1) , r.reactand_emus ) )
end
function Base.show(io::IO,r::SpeciesEMURxn)
    print(io,"$(r.flux): ")
    for si in r.reactand_emus
        print(io, "$(si) " )
    end
    print(io," --> ")
    print(io,r.product_emu)
end


# compares the flux and the product emu, this should uniquely identify the emurxn
function Base.isequal(r1::EMURxn,r2::EMURxn)
    return r1.flux==r2.flux && r1.product_emu==r2.product_emu
end
function Base.hash(r::EMURxn)
    return hash( string(r.flux, "$(hash(r.product_emu))" ) )
end

"""
compute_emu_decomp(n::Network,targets::Vector{Set{SAS}})
returns Vector{EMURxn}

computes EMU decomposition for the given network and the given target SpeciesEMUs
"""
function compute_emu_decomp(n::Network,targets::Vector{Set{SAS}})

    debuglevel=1

    # prepare:
    # compute mappings SAS <-> RAS , SAS <-> PAS
    all_sas::Set{SAS} = sas(n)
    all_ras::Set{RAS} = ras(n)
    all_pas::Set{PAS} = pas(n)
    d_sas2pas::Dict{SAS,Set{PAS}} = Dict{SAS,Set{PAS}}()
    d_pas2sas::Dict{PAS,SAS}      = Dict{PAS,SAS}()
    d_sas2ras::Dict{SAS,Set{RAS}} = Dict{SAS,Set{RAS}}()
    d_ras2sas::Dict{RAS,SAS}      = Dict{RAS,SAS}()
    for si in all_sas
        m_ras::Set{RAS} = filter( x1 -> x1.sas==si , all_ras )
        d_sas2ras[si] = Set(m_ras)
        for ri in m_ras
            d_ras2sas[ri] = si
        end
        m_pas::Set{PAS} = filter( x1 -> x1.sas==si , all_pas )
        d_sas2pas[si] = Set(m_pas)
        for pi in m_pas
            d_pas2sas[pi] = si
        end
    end

    if(debuglevel>=2)
        println("\n------------------------------------------------------\n")
        println(d_sas2pas)
        println("\n------------------------------------------------------\n")
        println(d_sas2ras)
        println("\n")
    end

    ## Start with decomposition:

    processed_emus::Set{SpeciesEMU} = Set{SpeciesEMU}()
    active_emus::Set{SpeciesEMU} = Set{SpeciesEMU}()

    ## check if we have target semus which span more than one metabolite, those will be split..
    for ti in targets
        if(!isinsinglespecies(ti))
            printf("Encountered multi-species target SEMUs, those will be split..")
            tis = splitSEMUIntoSpecies(target)
            for tisi in values(tis)
                push!( active_emus , tisi )
            end
        else
            push!( active_emus , ti )
        end
    end



    emu_rxns::Vector{EMURxn} = Vector{EMURxn}();

    while( ! isempty(active_emus) )
        active_emus = setdiff(active_emus,processed_emus)
        if(isempty(active_emus))
            break # we' 're done
        end

        # get next emu
        emu = pop!(active_emus)
        push!(processed_emus,emu)

        if(debuglevel>=2); println(string("Process EMU ",emu)); end

        # 1. construct EMU Reactions which gives rise to this SEMU
        # 1.1 find all PAS which contain this Species EMU
        pemus::Set{PAS} = findPASsForSEMU( d_sas2pas , emu )
        # 1.2 split into "products", product is (flux,species,ppos)
        sorted_pemus::Dict{Tuple{String,String,Int64},Set{PAS}} = splitPASsIntoProducts(pemus)

        # 2 loop over sorted_pemus and construct emu rxn for each

        for pemu_p::Tuple{String,String,Int64} in keys(sorted_pemus)
            if( isinflux( n.fluxes_db[pemu_p[1]] ) )
                if(debuglevel>=2); println("Hit Influx $(pemu_p[1])"); end
                continue
            end

            pemu::Set{PAS} = sorted_pemus[pemu_p]
            flux_i = n.fluxes_db[pemu_p[1]]
            rtransition = map(reverse, flux_i.transition)

            # 2.1 find all ras:
            # sort pemu!..
            vpemu::Vector{PAS} = sort(collect(pemu))
            vras::Vector{RAS} = map( x1 -> rtransition[x1] , vpemu )
            # 1.3 Split REMUs into reactands, "reactand" is (flux,species,rpos)
            sorted_ras::Dict{Tuple{String,String,Int64},Set{RAS}} = splitRASsIntoReactands(Set(vras))

            # create EMU Rxn..
            if(debuglevel>=2)
                println("Make EMU Rxn:")
                println("Target: \n $(vpemu)")
                println("REMUs Sizes: $(map( x -> length(sorted_ras[x]),keys(sorted_ras))) ")
                println("REMUs:")
                for ki in keys(sorted_ras)
                    println(sorted_ras[ki])
                end
                # \n  $(sorted_ras)\n")
            end

            # create rxn
            push!(emu_rxns, EMURxn( pemu_p[1] , pemu , collect(values(sorted_ras)) ) )

            # add mew semus to be processed:
            for ki in keys(sorted_ras)
                push!( active_emus , remu2semu(sorted_ras[ki]) )
            end
        end
    end
    return emu_rxns
end

#return Dict{String,Set{SAS}}
function splitSEMUIntoSpecies(sass::Set{SAS})
    all_species = map( x::SAS -> (x.species) , sass )
    species = Set(all_species)
    sortedspecies = Dict{String,Set{SAS}}()
    for si in species
        sortedspecies[si] = Set( filter( x1 -> x1.species==si , sass) )
    end
    return sortedspecies
end

function splitPASsIntoFluxes(sp::Set{PAS})
    pfluxes = Set(map( x::PAS -> x.flux , pemus ))
    sortedpas = Dict{String,Set{PAS}}()
    for pi in pfluxes
        sortedpas[si] = Set( filter( x1 -> x1.flux==pi , sp) )
    end
    return sortedpas
end

# product is ( "Flux" , "Species" , "ppos" )
function splitPASsIntoProducts(sp::Set{PAS})
    products::Set{Tuple{String,String,Int64}} =  Set{Tuple{String,String,Int64}}( map( x -> (x.flux,x.sas.species,x.ppos) , sp ) )
    sorted_pas = Dict{Tuple{String,String,Int64},Set{PAS}}()
    for pi in products
        sorted_pas[pi] = Set( filter( x1 -> ( (x1.flux==pi[1]) && (x1.sas.species==pi[2]) && (x1.ppos==pi[3]) ), sp) )
    end
    return sorted_pas
end



function isinsingleflux(sp::Set{RAS})
    if(isempty(sp)); return true; end
    sp1 = collect(sp)[1].flux
    return all( map( x -> x.flux == sp1 , sp ) )
end

# ! Function assumes that all are in same flux !
# reactand is ( "Flux" , "Species" , "rpos" )
function splitRASsIntoReactands(sp::Set{RAS})
    if(!isinsingleflux(sp)); error("::findProductEMUsForSEMU : predondition violated"); end
    reactands::Set{Tuple{String,String,Int64}} =  Set{Tuple{String,String,Int64}}( map( x -> (x.flux,x.sas.species,x.rpos) , sp ) )
    sorted_ras::Dict{Tuple{String,String,Int64},Set{RAS}} = Dict{Tuple{String,String,Int64},Set{RAS}}()
    for ri in reactands
        sorted_ras[ri] = Set( filter( x1 -> ( (x1.flux==ri[1]) && (x1.sas.species==ri[2]) && (x1.rpos==ri[3]) ), sp) )
    end
    return sorted_ras
end





# assumption: emu is inside one species!
function findPASsForSEMU( d_sas2pas::Dict{SAS,Set{PAS}} , emu::SpeciesEMU )
    # sanity check:
    if(!isinsinglespecies(emu)); error("::findProductEMUsForSEMU : predondition violated"); end

    # all PAS :
    all_pas::Set{PAS} = reduce( (a,b) -> union(a,b) , map( x -> d_sas2pas[x] , emu ) )

    return all_pas
end
