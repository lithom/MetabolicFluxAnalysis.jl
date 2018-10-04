using LightGraphs

import Base.isequal

# vertex, either SpeciesEMU or EMURxn
abstract type AbstractV
end

struct VSEMU <: AbstractV
    semu::SpeciesEMU
end
function Base.hash(x::VSEMU)
    return hash(string("abc",hash(x.semu)))
end

struct VRxn <: AbstractV
    r::EMURxn
end
function Base.hash(x::VRxn)
    return hash(string("cde",hash(x.r)))
end


function Base.isequal(x::AbstractV,y::AbstractV)
    #println("\n\nEquality Test for AbstractV")
    #println(typeof(x))
    #println(typeof(x))
    if( typeof(x) != typeof(y) )
        return false
    end
    if( typeof(x) == VSEMU )
        return isequal( x.semu , y.semu )
    end
    if( typeof(x) == VRxn )
        return isequal( x.r , y.r )
    end
end

mutable struct AVGraph
    g::SimpleDiGraph
    graph_db::Dict{Int64,AbstractV}
    graph_db_inv::Dict{AbstractV,Int64}
    last_index::Int64
end

function get_all_emurxns(g::AVGraph)
    vrxns = filter( x1 -> typeof(x1) == VRxn , collect( keys( g.graph_db_inv ) ) )
    rxns::Vector{EMURxn} = map( x1 -> x1.r , vrxns )
end

function get_all_semus(g::AVGraph)
    vsemus = filter( x1 -> typeof(x1) == VSEMU , collect( keys( g.graph_db_inv ) ) )
    semus::Vector{SpeciesEMU} = map( x1 -> x1.semu , vsemus )
end


function AVGraph()
    #return AVGraph(graph(Int64[],Edge{Int64}[]),Dict{Int64,AbstractV}(),Dict{AbstractV,Int64}(),0)
    return AVGraph(SimpleDiGraph{Int64}(),Dict{Int64,AbstractV}(),Dict{AbstractV,Int64}(),0)
end
function av_add_vertex(g::AVGraph,v::AbstractV)
    #println("n vertices: $(nv(g.g))")
    if(haskey(g.graph_db_inv,v))
        return
    else
        g.last_index = g.last_index+1
        add_vertex!( g.g )
        g.graph_db[g.last_index] = v
        g.graph_db_inv[v] = g.last_index
    end
end
function av_add_edge(g::AVGraph, v1::AbstractV , v2::AbstractV)
    av1 = g.graph_db_inv[v1]
    av2 = g.graph_db_inv[v2]
    add_edge!(g.g,av1,av2)
end
function split_into_connected_components(g::AVGraph)
    split_components = Vector{AVGraph}()
    wccs = weakly_connected_components( g.g )

    for wcci in wccs
        (g_i,vmap) = induced_subgraph(g.g,wcci)
        new_dict::Dict{Int64,AbstractV} = Dict{Int64,AbstractV}()
        new_dict_rev::Dict{AbstractV,Int64} = Dict{AbstractV,Int64}()
        for zi in 1:length(wcci)
            new_dict[ zi ] = g.graph_db[ vmap[zi] ]
            new_dict_rev[ g.graph_db[ vmap[zi] ] ] = zi
        end
        push!( split_components , AVGraph(g_i,new_dict,new_dict_rev,length(wcci)) )
    end

    return split_components
end


"""
create_emu_graph(rxns::Vector{EMURxn})
returns Dict{Int64,Vector{AVGraph}}

# what it does:
  1. it sorts emu rxns according to size of target EMUs
  2. it computes the connected compoments
  it then returns the connected components, sorted by size
"""
function create_emu_graph(rxns::Vector{EMURxn})

    sizes::Set{Int64} = Set{Int64}( map( x1 -> length(x1.product_emu) , collect(rxns) ) )
    sorted_rxns      = Dict{Int64,Set{EMURxn}}()

    network_components::Dict{Integer,Vector{AVGraph}} = Dict{Integer,Vector{AVGraph}}()

    for si = sizes
        rxi = filter( x1 -> (length(x1.product_emu)==si) , rxns )
        sorted_rxns[si] = Set(rxi)

        # create graphs and find connected compoments..
        #g::GenericGraph{AbstractV,Edge{AbstractV}} = graph(AbstractV[],Edge{AbstractV}[])
        #g::AVGraph(simple_graph(0),Dict{Int64,AbstractV}(),Dict{AbstractV;Int64}(),0)
        g::AVGraph = AVGraph()

        # add all rxns
        for ri in rxi
                semu_p  = pemu2semu( ri.product_emu )
                semus_r = map( xi -> remu2semu(xi) , ri.reactand_emus )

                # add vertices and edges
                vp = VSEMU(semu_p)
                # add_vertex!( g , vp )
                av_add_vertex(g,vp)
                vr = VRxn(ri)
                #add_vertex!( g , vr )
                av_add_vertex(g,vr)
                #add_edge!(g,vr,vp)
                av_add_edge(g,vr,vp)

                v_sr = map( xi -> VSEMU(xi) , semus_r )
                #for vsri in v_sr ; add_vertex!( g , vsri ) ; end
                #for vsri in v_sr ; add_edge!( g , vsri , vr ) ; end
                for vsri in v_sr ; av_add_vertex( g , vsri ) ; end
                for vsri in v_sr ; av_add_edge( g , vsri , vr ) ; end
        end

        #print(g.g)
        #return g

        # now split into connected compoments
        g_split::Vector{AVGraph} = split_into_connected_components(g)
        network_components[si] = g_split
    end
    return network_components
end

"""
get_roots_of_av_graph(g::AVGraph)
returns Tuple{Vector{SpeciesEMU},Vector{SpeciesEMU}}

returns all semu vertices with zero indegree, and the complement
"""
function get_roots_of_av_graph(g::AVGraph)
    roots::Vector{SpeciesEMU} = Vector{AbstractV}()
    non_roots::Vector{SpeciesEMU} = Vector{AbstractV}()
    for vi::AbstractV in keys(g.graph_db_inv)
        if( vi isa VSEMU )
            if( length( inneighbors( g.g , g.graph_db_inv[vi] ) )==0)
                push!(roots,vi.semu)
            else
                push!(non_roots,vi.semu)
            end
        end
    end
    return ( roots , non_roots )
end
