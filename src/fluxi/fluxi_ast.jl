
import Base.isequal
import Base.hash
import Base.show

"""
abstract type ASTNumNode
    represents the AST nodes which evaluate to a number, when supplied with
    all the data (mds, fluxes)
"""
abstract type ASTNumNode end

struct ASTMatrix
    m::Array{ASTNumNode,2}
end

"""
ASTSymmOp
has fields: op::Char , args::Vector{ASTNumNode}
represents addition and multiplication operations
"""
struct ASTSymmOp <: ASTNumNode
    op::Char
    args::Vector{ASTNumNode}
end

"""
ast_add(a::ASTNumNode,b::ASTNumNode)
creates sum of two ASTNumNode objects
"""
function ast_add(a::ASTNumNode,b::ASTNumNode)
    return ASTSymmOp('+',[a,b])
end

function ast_add(a::ASTSymmOp,b::ASTNumNode)
    if(a.op=='+')
        push!( a.args , b )
        return a
    else
        return ASTSymmOp('+',[a,b])
    end
end
function ast_add(a::ASTNumNode,b::ASTSymmOp)
    ast_add(b,a)
end

"""
ASTRealValue has field val::Float64
and represents a hard-coded real value
"""
struct ASTRealValue <: ASTNumNode
    val::Float64
end

"""
ASTEMUMassRatio has fields
  semu::SpeciesEMU , mass:Float64

and represents a mass ratio (which has to be supplied to compute the FW-Sim)
"""
struct ASTEMUMassRatio <: ASTNumNode
    semu::SpeciesEMU
    mass::Float64 # here we could also have int, basically..
end
function Base.hash(emr::ASTEMUMassRatio)
    mass = round(emr.mass)
    return xor( Base.hash(emr.semu) , Base.hash(mass) )
end
function Base.isequal(a::ASTEMUMassRatio,b::ASTEMUMassRatio)
    return Base.isequal(a.semu,b.semu) && isequal(round(a.mass),round(b.mass))
end
function Base.isless(a::ASTEMUMassRatio,b::ASTEMUMassRatio)
    if( isless( a.semu , b.semu ) )
        return true
    end
    if( isequal(a.semu,b.semu) )
        return round(a.mass) < round( b.mass )
    end
    return false
end
function Base.show(io::IO,emr::ASTEMUMassRatio)
    print(io,"ASTEMR:$(emr.semu):M=$(round(emr.mass))")
end

"""
ASTFluxRate has fields
  flux::String

and represents a flux rate (which has to be supplied to compute the FW-Sim)
"""
struct ASTFluxRate <: ASTNumNode
    flux::String
end
function Base.hash(f::ASTFluxRate)
    return Base.hash("flux$(f.flux)")
end
function Base.isequal(a::ASTFluxRate,b::ASTFluxRate)
    return isequal(a.flux,b.flux)
end
function Base.isless(a::ASTFluxRate,b::ASTFluxRate)
    return a.flux < b.flux
end


"""
ASTEMUSimStep has fields
  A::ASTMatrix , X::ASSTMatrix , B::ASTMatrix , Y::ASSTMatrix

and represents a step of the FW-Sim. Notation is according to Antoniewicz2006,
i.e. new SEMU masses are in matrix X, and given ones are in matrix Y,
i.e we have to compute the step X = inv(A) * B * Y
"""
struct ASTEMUSimStep
    A::ASTMatrix
    X::ASTMatrix
    B::ASTMatrix
    Y::ASTMatrix
end

"""
ASTEMUSim has fields
  sim_steps::Dict{Integer,Vector{ASTEMUSimStep}}
  configured_semus::Vector{SpeciesEMU}

  and represents all information required to compute the FW-Sim, for given fluxes
  and mds for the semus in .configured_semus
  the simulation steps are sorted by size of EMUs
"""
struct ASTEMUSim
    sim_steps::Dict{Integer,Vector{ASTEMUSimStep}}
    configured_semus::Vector{SpeciesEMU}
end


"""
create_ast_from_emu_graph( eg::Dict{Integer,Vector{AVGraph}} )
returns ASTEMUSim

computes the AST which represents the forward simulation
"""
function create_ast_from_emu_graph( eg::Dict{Integer,Vector{AVGraph}} )
    required_semus::Set{SpeciesEMU} = Set{SpeciesEMU}()
    determined_semus::Set{SpeciesEMU} = Set{SpeciesEMU}()
    for si in sort( collect( keys(eg) ) ) # we have to go from small to large, as in the emu decomp..
        for gi in eg[si]
            ( roots::Vector{SpeciesEMU} , non_roots::Vector{SpeciesEMU} ) = get_roots_of_av_graph(gi)
            req_new = setdiff( Set( roots ) , determined_semus )
            union!( required_semus , req_new )
            union!( determined_semus , Set(non_roots ) )
        end
    end
    print("\nidentified $(length(required_semus)) required semus:\n")
    for rsi in required_semus; print("  $(rsi)\n"); end
    # ok, unroll them into single sas..
    required_sas = reduce( (x,y) -> union(x,y) , required_semus )
    #required_sas = sort(required_sas)
    print("\nSAS Set: $(length(required_sas)) required sas:\n")
    for sasi in sort(collect(required_sas))
        print(" $(sasi) ")
    end
    print("\n")

    return create_ast_from_emu_graph( eg , required_sas)
end


# configured_sas : all SAS which are given by the substrate configuration
function create_ast_from_emu_graph( eg::Dict{Integer,Vector{AVGraph}} , configured_sas::Set{SAS})

    debug_level = 1;

    sizes = sort( collect( keys( eg ) ) )

    # for sanity we keep track of all determined semu on the way..
    all_determined_semus::Set{SpeciesEMU} = Set{SpeciesEMU}()

    ast_emu_sim = ASTEMUSim(Dict{Integer,Vector{ASTEMUSimStep}}() , Vector{SpeciesEMU}())

    for si in sizes
        ast_emu_sim.sim_steps[si] = Vector{ASTEMUSimStep}()
        ni_cnt = 0
        for ni in eg[si]
            ni_cnt += 1
            # ok, lets build the matrix..
            # Step 1: partition semus into determined / target
            # Explanation: determined are:
            #    1. all non-max-size EMU ,
            #    2. EMUs containing only configured SAS
            # OR: easier: undetermined EMUs are all the product EMUs
            s_det::Set{SpeciesEMU}    = Set{SpeciesEMU}()
            s_nondet::Set{SpeciesEMU} = Set{SpeciesEMU}()

            emu_rxns  = get_all_emurxns(ni)
            semu_rxns = map( x1 -> SpeciesEMURxn(x1) , emu_rxns )
            semus     = get_all_semus(ni)

            for sri in semus
                if( length(sri) < si )
                    push!( s_det , sri )
                else
                    if( issubset( sri , configured_sas ) )
                        push!( s_det , sri )
                    else
                        push!( s_nondet , sri )
                    end
                end
            end

            if(debug_level>=1)
                print(  "--------------------------------------------------------------\n")
                print(  "NETWORK OF SIZE: $(si) , $(ni_cnt)/$(length(eg[si])) \n")
                print(  "--------------------------------------------------------------\n")
                print(  "ALL EMU RXNs: ------------------------------------------------\n")
                for rxi in semu_rxns; print("  $(rxi)\n"); end
                print(  "SEMUs (determined): ------------------------------------------\n")
                #for si in semus; print("  $(si)\n"); end
                for sid in s_det; print("  $(sid)\n"); end
                print(  "SEMUs (non-determined): --------------------------------------\n")
                #for si in semus; print("  $(si)\n"); end
                for sind in s_nondet; print("  $(sind)\n"); end
                print(  "--------------------------------------------------------------\n\n")
            end

            ## NOW: Build the matrices:
            #  X : Rows correspond to the non-det EMUs
            #      Cols correspond to the weight
            #  Y : Rows correspond to the determined EMUs
            #      Cols correspond to the weight
            #  A : Rows correspond to the non-det EMUs
            #      Cols correspond to the non-det EMUs
            #  B : Rows correspond to the non-det EMUs
            #      Cols correspond to the det EMUs

            # Order of rows / columns determined by:
            nondet_semu_order = sort(collect(s_nondet))
            det_semu_order    = sort(collect(s_det))

            ndspos = (x1::SpeciesEMU) -> find( xi -> xi==x1 , nondet_semu_order )
            #dspos  = (x1::SpeciesEMU) -> find( xi -> xi==x1 , det_semu_order )

            # HOW TO BUILD:
            # For each EMU we consider the influxes.
            # The influxes times the EMU mixture has to be equal
            # to the combined influxes of the EMU mixtures which feed the EMU

            # Build A, that's easy;
            A::Array{ASTNumNode,2} = Array{ASTNumNode,2}(length(nondet_semu_order),length(nondet_semu_order))
            for ns in nondet_semu_order
                #build this row:

                # find all emu rxns with product ns
                emu_rxns_ns = filter( er -> isequal( pemu2semu(er.product_emu) , ns ) , emu_rxns )
                # create the "entry on the diagonal"
                diag_entry::ASTNumNode = ASTSymmOp('+',Vector{ASTNumNode}())
                for eri in emu_rxns_ns
                    push!( diag_entry.args , ASTSymmOp('*',[ASTRealValue(-1),ASTFluxRate(eri.flux)] ) )
                end
                A[ndspos(ns),ndspos(ns)] = diag_entry
                # create the other entries:
                for eri in emu_rxns_ns
                    # skip for emu_rxns with more than one remu (becuase those cannot have undetermined remus..)
                    if( length( eri.reactand_emus ) > 1 ); continue; end
                    rpos  = ndspos( remu2semu( eri.reactand_emus[1] ) )
                    if(!isempty(rpos))
                        ndpos = ndspos(ns) ; @assert length(ndpos)==1 ; ndpos = ndpos[1] ;
                        @assert length(rpos)==1 ; rpos = rpos[1] ;
                        if( ! isassigned(A,ndpos,rpos))
                            A[ndpos,rpos] = ASTFluxRate(eri.flux)
                        else
                            A[ndpos,rpos] = ast_add( A[ndpos,rpos] , ASTFluxRate(eri.flux) )
                        end
                    end
                end
            end
            if(debug_level>=2) ; print(A) ; end

            # Build X (trivial)
            #
            if(debug_level>=2); println("\n\nsi : $(si)\n\n") ; end
            X::Array{ASTNumNode,2} = Array{ASTNumNode,2}(length(nondet_semu_order),si+1)
            for zxi in 0:si
                for ns in nondet_semu_order
                    X[ ndspos(ns) ,zxi+1] = ASTEMUMassRatio(ns,zxi)
                end
            end


            # Build B, similar to A, just with cols are the emu rxns with det.
            # reactand emus. NOTE: if one reactand is det. then all are det. ( as all are non-maximal size)

            det_rxns_1 = filter( x1 -> isempty(setdiff( x1.reactand_emus , s_det )) , semu_rxns )
            det_rxns_2 = filter( x1 -> !isempty(intersect( x1.reactand_emus , s_det )) , semu_rxns )
            @assert length(det_rxns_1)==length(det_rxns_2)

            det_rxns      = det_rxns_1
            det_rxn_order = sort(det_rxns) # this order defines the rows of B and the cols of Y
            drpos         = (x1::SpeciesEMURxn) -> find( xi->xi==x1 , det_rxn_order )[1]

            B::Array{ASTNumNode,2} = Array{ASTNumNode,2}(length(nondet_semu_order),length(det_rxns))
            Y::Array{ASTNumNode,2} = Array{ASTNumNode,2}(length(det_rxns),si+1)
            # Build B
            for ns in nondet_semu_order
                for rd in det_rxn_order
                    if( isequal( rd.product_emu , ns ) )
                        # !! MINUS
                        B[ndspos(ns),drpos(rd)] = ASTSymmOp('*',[ASTRealValue(-1),ASTFluxRate(rd.flux)] ) #ASTFluxRate(rd.flux)
                    end
                end
            end
            # Build Y , colums are masses 0:si
            for rd in det_rxn_order
                for mi in 0:si
                    # now collect all combinations of masses to get mi
                    rd_semus = rd.reactand_emus
                    rd_semus_sizes = map( rxi -> length(rxi) ,  rd_semus)
                    mmcs = matching_mass_combinations( mi , rd_semus_sizes)
                    yi::ASTNumNode = ASTSymmOp('+',Vector{ASTNumNode}(0))
                    for mmcsi in mmcs
                        yi_p::ASTNumNode = ASTSymmOp('*',Vector{ASTNumNode}(0))
                        wmicnt = 0
                        for wmi in mmcsi
                            wmicnt += 1
                            push!( yi_p.args , ASTEMUMassRatio( rd_semus[wmicnt] , wmi ) )
                        end
                        push!(yi.args,yi_p)
                    end
                    Y[drpos(rd),mi+1] = yi
                end
            end

            emusimstep_i = ASTEMUSimStep(ASTMatrix(A),ASTMatrix(X),ASTMatrix(B),ASTMatrix(Y))
            push!( ast_emu_sim.sim_steps[si] , emusimstep_i )
        end
    end

    if(debug_level>=1) ; println("\ndone..\n") ; end
    return ast_emu_sim
end


function construct_all_mass_ratios_for_measurement(fd::FluxiData , measurement::String )
    construct_all_mass_ratios_for_measurement(fd.measurements[measurement])
end

function construct_all_mass_ratios_for_measurement(vmd::Vector{MassDistributionMeasurement})
    smrs::Set{ASTEMUMassRatio} = Set{ASTEMUMassRatio}()
    for mdi in vmd
        #print(mdi)
        semu = mdi.sass
        for mi=0:length(semu)
            push!(smrs,ASTEMUMassRatio(semu,mi))
        end
    end
    return sort(collect(smrs))
end

function construct_all_mass_ratios_for_measurement(semus::Vector{SpeciesEMU})
    smrs::Set{ASTEMUMassRatio} = Set{ASTEMUMassRatio}()
    for si in semus
        for mi=0:length(si)
            push!(smrs,ASTEMUMassRatio(si,mi))
        end
    end
    return sort(collect(smrs))
end

function construct_mds(fd::FluxiData,substrate_conf::String,mds::Vector{ASTEMUMassRatio})
    construct_mds( fd.substrate_conf[substrate_conf] , mds )
end


"""
construct_mds(vmd::Vector{ConfiguredSEMUMixture},mds::Vector{ASTEMUMassRatio})
returns Vector{Float64}
computes the mds which are required to run the fw-sim from a ConfiguredSEMUMixture
object.

# Note: Required to init the fw simulation
"""
function construct_mds(vmd::Vector{ConfiguredSEMUMixture},mds::Vector{ASTEMUMassRatio})
    md_values = Vector{Float64}(length(mds))
    # preprocess the vmd..
    d_cs = Dict{ConfiguredSEMU,Float64}()
    for vmdi in vmd; d_cs = merge(d_cs,vmdi); end

    cnt_md = 0
    for md in mds
        cnt_md +=1
        @assert ( isinsinglespecies( md.semu ) )
        atom_positions = sort( map( s1 -> s1.pos , collect(md.semu) ) )

        # find the entries in d_cs which cover this semu..
        csemus = filter( x1 -> isempty( setdiff( md.semu , x1.semu ) ) , collect(keys(d_cs)) )
        #println("\nTotal values: $(sum( map( x1 -> d_cs[x1] , csemus ) ))\n")
        @assert sum( map( x1 -> d_cs[x1] , csemus ) ) > 0.999
        @assert sum( map( x1 -> d_cs[x1] , csemus ) ) < 1.001

        # now just select all which have md.mass labeled atoms and sum up..
        m_semus = filter( x1 -> round( sum(x1.config[atom_positions]))==round(md.mass) , csemus )
        if(isempty(m_semus))
            md_values[cnt_md] = 0.0;
        else
            md_values[cnt_md] = sum( map( x -> d_cs[x] , m_semus ) )
        end
    end
    return md_values
end


function get_mass_distribution_measurement_value( fd::FluxiData , mname::String , emu::SpeciesEMU , mass::Float64 )
    return get_mass_distribution_measurement_value( fd , mname , ASTEMUMassRatio(emu,mass) )
end

function get_mass_distribution_measurement_value( fd::FluxiData , mname::String , edr::ASTEMUMassRatio )
    mds::Array{MassDistributionMeasurement,1} = fd.measurements[mname]
    # 1. find the right measurement:
    mdsi = filter( x -> x.sass == edr.semu , mds)
    if(length(mdsi)==0);
        println("\nNo measurement found for semu $(edr.semu)")
        return NaN;
    end
    if(length(mdsi)>1)
        println("\n!!MULTIPLE MEASUREMENTS FOUND!! for semu $(edr.semu)")
        data = Array{Float64,2}(3,length(mdsi))
        for zi in 1:length(mdsi)
            data[:,zi] = mdsi[zi].data[edr.mass]
        end
    else
        data = Array{Float64,2}(3,1)
        data[:,1] = mdsi[1].data[edr.mass]
    end
    return data
end




# AST iterator:
function ast_traverse_df!( n::ASTEMUSim , v::Vector{Any})
    for si in sort(collect(keys(n.sim_steps)))
        for ssi in n.sim_steps[si]
            ast_traverse_df!( ssi.A , v )
            ast_traverse_df!( ssi.X , v )
            ast_traverse_df!( ssi.B , v )
            ast_traverse_df!( ssi.Y , v )
            push!( v , ssi )
        end
    end
    push!( v , n )
end
function ast_traverse_df!( n::ASTMatrix , v::Vector{Any})
    for zi=1:size(n.m,2)
        for zj=1:size(n.m,1)
            if(isassigned(n.m,zj,zi))
                ast_traverse_df!( n.m[zj,zi] , v )
            end
        end
    end
    push!( v , n )
end
function ast_traverse_df!( n::ASTSymmOp , v::Vector{Any})
    for ai in n.args
        ast_traverse_df!( ai , v )
    end
    push!( v , n )
end
function ast_traverse_df!( n::ASTNumNode , v::Vector{Any})
    push!( v , n )
end
