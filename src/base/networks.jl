
import Base.isless

"""
Species has fields id::String and na::Int64
and represents a metabolite/species
"""
struct Species
    id::String
    na::Int64 # number of atoms
end

function Species(id::String)
    return Species(id,0)
end
function isless(x::Species,y::Species)
    return x.id < y.id
end

"""
Flux has fields:
  id::String , E::LC{String} , P::LC{String} , transition::Dict{RAS,PAS}
and represents a flux/reaction
"""
struct Flux
    id::String
    E::LC{String}
    P::LC{String}
    transition::Dict{RAS,PAS}

    function Flux(id::String,E::LC{String},P::LC{String},transition::Dict{RAS,PAS})
        E2 = E
        P2 = P
        if( any( collect(values(E)).<0) || any( collect(values(P)).<0) )
            error("Flux Constructor : Stoich. with neg. entries..")
        end
        if( isempty(E) && isempty(P) )
            error("Flux Constructor : Stoich. without any entries..")
        end
        if( isempty(E) )
            E2 = LC{String}( "*" => 1.0)
        end
        if( isempty(P) )
            P2 = LC{String}( "*" => 1.0)
        end
        new(id,E2,P2,transition)
    end
end


function Flux(id::String,stoichiometry::LC{String})
    # auto split stoich. into positive and negtive parts:
    lce = LC{String}()
    lcp = LC{String}()
    for ki in keys(stoichiometry)
        if(stoichiometry[ki]>0)
            lcp[ki] = stoichiometry[ki]
        else
            lce[ki] = -stoichiometry[ki]
        end
    end
    return Flux(id,lce,lcp,Dict{RAS,PAS}())
end

function isless(x::Flux,y::Flux)
    return x.id < y.id
end
function isioflux(f::Flux)
    return ( isinflux(f) || isoutflux(f) )
end
function isinflux(f::Flux)
    return haskey(f.E, "*")
end
function isoutflux(f::Flux)
    return haskey(f.P, "*")
end

function create_reverse_flux( f::Flux )
    name_bw_flux = string(f.id,"_bw")
    return create_reverse_flux( f , name_bw_flux )
end
function create_reverse_flux( f::Flux , name_bw_flux::String)
    tmap_reverse = Dict{RAS,PAS}()
    for ri in keys(f.transition)
        pi = f.transition[ri]
        tmap_reverse[ RAS(name_bw_flux , pi.sas , pi.ppos) ] = PAS(name_bw_flux,ri.sas,ri.rpos)
    end
    return Flux(name_bw_flux,f.P,f.E,tmap_reverse)
end


"""
Network has fields:
  fluxes::Vector{String} , species::Vector{String}
  fluxes_db::Dict{String,Flux} , species_db::Dict{String,Species}
and represents a flux/reaction
"""
struct Network
    fluxes::Vector{String}
    species::Vector{String}
    fluxes_db::Dict{String,Flux}
    species_db::Dict{String,Species}
end


function Network()
    Network( Vector{Flux}() )
end

"""
sas(id_metabolite::String,n::Network)
returns Vector{SAS}
computes all SAS for a given metabolite
"""
function sas(ss::String,n::Network)
    vsas::Vector{SAS} = Vector{SAS}()
    s::Species = n.species_db[ss]
    print(s)
    for i in 1:s.na
        push!( vsas , SAS(ss,i) )
    end
    return vsas
end

"""
ras(id_flux::String,n::Network)
returns Vector{RAS}
computes all RAS for a given flux
"""
function ras(fs::String,n::Network)
    f::Flux = n.fluxes_db[fs]

    if(isinflux(f)); return Vector{RAS}(); end

    e::Vector{String} = collect( keys( f.E ) )
    sort(e)
    vras::Vector{RAS} = Vector{RAS}()
    for ei::String in e
        if(n.species_db[ei].na>0)
            nrf  = f.E[ei]
            nr::Int64 = Int64(round(nrf))
            if( abs(nrf-nr) > (1e-8) )
                error("Isotopomer Species with non-integer number of reactands..")
            end
            for ri in 1:nr
                for rsas in sas(ei,n)
                    push!( vras , RAS( fs , rsas , ri ) )
                end
            end
        end
    end
    return vras
end

"""
pas(id_flux::String,n::Network)
returns Vector{PAS}
computes all PAS for a given flux
"""
function pas(fs::String,n::Network)
    f::Flux = n.fluxes_db[fs]

    if(isoutflux(f)); return Vector{PAS}(); end

    p::Vector{String} = collect( keys( f.P ) )
    sort(p)
    vpas::Vector{PAS} = Vector{PAS}()
    for pi::String in p
        if(n.species_db[pi].na>0)
            nrf  = f.P[pi]
            nr::Int64 = Int64(round(nrf))
            if( abs(nrf-nr) > (1e-8) )
                error("Isotopomer Species with non-integer number of reactands..")
            end
            for ri in 1:nr
                for rsas in sas(pi,n)
                    push!( vpas , PAS( fs , rsas , ri ) )
                end
            end
        end
    end
    return vpas
end

"""
ras(n::Network)
returns Vector{SAS}
computes all SAS of all species of the network
"""
function sas(n::Network)
    all_sas::Set{SAS} = Set{SAS}()
    for si in n.species
        all_sas = union( all_sas , Set(sas(si,n)) )
    end
    return all_sas
end

function sas(csemu::ConfiguredSEMU)
    return Set( csemu.semu )
end



"""
ras(n::Network)
returns Vector{RAS}
computes all RAS of all fluxes of the network
"""
function ras(n::Network)
    all_ras::Set{RAS} = Set{RAS}()
    for fi in n.fluxes
        all_ras = union( all_ras , Set(ras(fi,n)) )
    end
    return all_ras
end

"""
pas(n::Network)
returns Vector{PAS}
computes all PAS of all fluxes of the network
"""
function pas(n::Network)
    all_pas::Set{PAS} = Set{PAS}()
    for fi in n.fluxes
        all_pas = union( all_pas , Set(pas(fi,n)) )
    end
    return all_pas
end

"""
  get_N( n::Network ; sparse=false)
returns Array{Float64,2}
computes the stoichiometric matrix of the network
"""
function get_N(net::Network ; sparse_n=false)
    (n1,n2,n3) = get_full_N(net ; sparse_n=sparse_n)
    return n1
end

"""
  get_full_N( n::Network ; sparse=false)
returns Tuple{ Array{Float64,2} , Array{Float64,2}  , Array{Float64,2} }
computes the stoichiometric matrix of the network result[1], and the positive / negative
entries of the stoichiometric matrix (result[2]/result[3])
"""
function get_full_N(net::Network ; sparse_n=false)

    n_plus::Array{Float64,2}  = zeros(length(net.species),length(net.fluxes))
    n_minus::Array{Float64,2} = zeros(length(net.species),length(net.fluxes))

    d_s::Dict{String,Int64}   = Dict{String,Int64}()
    for zi=1:length(net.species); d_s[net.species[zi]] = zi; end

    for zi=1:length(net.fluxes)
        fi = net.fluxes_db[ net.fluxes[zi] ]
        for ki in keys(fi.E); if(ki=="*"); continue; end ; n_minus[d_s[ki],zi] = fi.E[ki]; end
        for ki in keys(fi.P); if(ki=="*"); continue; end ; n_plus[d_s[ki],zi]  = fi.P[ki]; end
    end

    n = (n_plus - n_minus)

    if(!sparse_n)
        return n,n_plus,n_minus
    else
        n_sp = SparseArrays.sparse(n); n_sp_p = SparseArrays.sparse(n_plus); n_sp_m = SparseArrays.sparse(n_minus);
        return n_sp , n_sp_p , n_sp_m
    end
end

"""
is_sx_flux(f::Flux)
returns Bool

checks if flux is sx flux, which is by convention equivalent of having the
"*" id as a key in reactand lc or product lc
"""
function is_sx_flux(f::Flux)
    return haskey( f.E , "*" ) || haskey( f.P , "*" )
end

"""
get_sx_fluxes(net::Network)
returns Vector{Flux}
returns sx fluxes of the network
"""
function get_sx_fluxes(net::Network)
    return map( x -> net.fluxes[x] , get_idx_sx_fluxes(net) )
end
function get_idx_sx_fluxes(net::Network)
    fsx = Vector{Int64}()
    fi_cnt = 0
    for fi in net.fluxes
        fi_cnt += 1
        if( is_sx_flux( net.fluxes_db[fi] ) ) ; push!(fsx,fi_cnt) ; end
    end
    return fsx
end
function get_non_sx_fluxes(net::Network)
    return map( x -> net.fluxes[x] , get_idx_non_sx_fluxes(net) )
end
function get_idx_non_sx_fluxes(net::Network)
    fsx = Vector{Int64}()
    fi_cnt = 0
    for fi in net.fluxes
        fi_cnt += 1
        if( !is_sx_flux( net.fluxes_db[fi] ) ) ; push!(fsx,fi_cnt) ; end
    end
    return fsx
end



"""
get_flux_polyhedron(net::Network , constraints::Vector{LinConstr} ; max_val=Inf)
returns Polyhedron

creates the flux polyhedron
CONVENTION:
System exchange fluxes are reversible, all other fluxes have 0 as lower bound
"""
function get_flux_polyhedron(net::Network , constraints::Vector{LinConstr} ; max_val=Inf)
     p_0 = get_flux_polyhedron_vanilla(net::Network , constraints::Vector{LinConstr} )

     d_s::Dict{String,Int64}   = Dict{String,Int64}()
     for zi=1:length(net.fluxes); d_s[net.fluxes[zi]] = zi; end

     # identify non-system-exchange fluxes
     fnsx = get_idx_non_sx_fluxes(net)
     #print(fnsx)

     lb =  Inf * -ones(length(net.fluxes))
     ub =  Inf *  ones(length(net.fluxes))
     lb[fnsx] = 0

     # if max_val, create bouding box. [0,max_val] for nsx fluxes, AND [-max_val,max_val] for sx fluxes
     if(max_val<Inf)
         ub[ub.>max_val] = max_val
         lb[ get_idx_sx_fluxes(net) ] = -max_val
     end

     return Polyhedron( p_0 , bounding_box(lb,ub) )
end


"""
get_flux_polyhedron_vanilla(net::Network , constraints::Vector{LinConstr} )
returns Polyhedron

creates the flux polyhedron, without >0 constraint on non-sx-fluxes

This function does not add the greater than zero constraints to the non-system-
exchange fluxes
"""
function get_flux_polyhedron_vanilla(net::Network , constraints::Vector{LinConstr} )
    m_n = get_N(net)
    p_c = lc2polyhedron(net , constraints)
    return Polyhedron( p_c.G , p_c.h , [m_n;p_c.A] , [zeros(size(m_n,1),1);p_c.b]  )
end

"""
lc2polyhedron(net::Network , lcs::Vector{LinConstr})
returns Polyhedron
creates polyhedron from the given constraints
"""
function lc2polyhedron(net::Network , lcs::Vector{LinConstr} )
    n = length(net.fluxes)
    d_s::Dict{String,Int64}   = Dict{String,Int64}()
    for zi=1:length(net.fluxes); d_s[net.fluxes[zi]] = zi; end
    d_s["1"] = n+1 # to collect the "1" data

    G::Array{ Float64 , 2} = Array{ Float64 }(0,n)
    h::Array{ Float64 , 2} = Array{ Float64 }(0,1)
    A::Array{ Float64 , 2} = Array{ Float64 }(0,n)
    b::Array{ Float64 , 2} = Array{ Float64 }(0,1)

    v_G =  Array{ RowVector{Float64,Array{Float64,1}},1 }()
    v_A =  Array{ RowVector{Float64,Array{Float64,1}},1 }()
    v_h =  Array{ Float64 , 1 }()
    v_b =  Array{ Float64 , 1 }()

    G_s = spzeros(0,n)
    A_s = spzeros(0,n)

    for lci in lcs
        row_lhs = resolve_lc( n+1 , d_s , lci.lhs )
        row_lhs_one = row_lhs[n+1]
        row_lhs = row_lhs[1:n]
        row_rhs = resolve_lc( n+1 , d_s , lci.rhs )
        row_rhs_one = row_rhs[n+1]
        row_rhs = row_rhs[1:n]

        # copy all to rhs, except for "one" part:
        p_row_lhs = row_lhs - 1.0 * row_rhs
        p_row_rhs = row_rhs_one - row_lhs_one

        if(lci.constr=='=')
            #A = [A;p_row_lhs']; b = [b;p_row_rhs]; continue;
            #push!( v_A , p_row_lhs'); push!( v_b, p_row_rhs); continue;
            A_s = vcat( A_s , p_row_lhs' ); push!( v_b, p_row_rhs); continue;
        end
        if(lci.constr=='<')
            #G = [G;p_row_lhs']; h = [h;p_row_rhs]; continue;
            #push!( v_G , p_row_lhs'); push!(v_h,p_row_rhs); continue;
            G_s = vcat( G_s , p_row_lhs' ); push!(v_h,p_row_rhs); continue;
        end
        if(lci.constr=='>')
            #G = [G;-p_row_lhs']; h = [h;-p_row_rhs]; continue;
            #vG = push!(-p_row_lhs'); push!(v_h,-p_row_rhs); continue;
            G_s = vcat( G_s , -p_row_lhs' ); push!(v_h,-p_row_rhs); continue;
        end
    end
    #print("\nDONE!\n")
    #println("\n\n $(G) \n\n $(h) \n\n $(A) \n\n $(b) \n\n ")
    #G = reduce( (x,y) -> [x;y] , v_G )
    #A = reduce( (x,y) -> [x;y] , v_A )
    h = [ h ; v_h]
    b = [ b ; v_b]
    G = full(G_s)
    A = full(A_s)
    return Polyhedron(G,h,A,b)
end
