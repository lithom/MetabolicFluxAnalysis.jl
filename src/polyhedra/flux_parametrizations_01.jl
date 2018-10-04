
using MetabolicFluxAnalysis

"""
  WX_Transform
represents a parametrization in terms of net / xch flux for reversible fluxes.

Usage to resolve indeces of the parametrization:
  get_wxflux_net(net,wxt,fid) : returns the position of the net flux for the
                                WX-pair which contains flux fid
  get_wxflux_xch(net,wxt,fid) : returns the position of the xch flux for the
                                WX-pair which contains flux fid

#NOTE: The indeces of the WX parametrization work as follows:
first, the wx-pair is resolved for a given flux, then the position of the FIRST
flux of that pair is always used to parametrize the net flux, the SECOND one
always to parametrize the XCH flux.
"""
struct WX_Transform
    wparams::Array{Int64,2} # wparams[i ,1/2] is index of ith net / sx flux (parametrization)
                             # by convention, the net flux is at the pos of the first pos. original flux,
                             # and the sx flux is at the second pos. original flux

    f::Any
end

function eval_wx_transform(wt::WX_Transform,xp::Array{Float64,2})
    xf = copy(xp)
    for zi in 1:size(wt.wparams,1)
        wti = wt.wparams[zi,:]
        for zx in 1:size(xp,2)
            (wfa::Float64,wfb::Float64) = wt.f( xp[ wti[1] , zx ] , xp[ wti[2] , zx ] )
            xf[ wti , zx ] = [wfa;wfb]
        end
    end
    return xf
end

"""
  get_wxflux_net(net::Network , wt::WX_Transform , fid::String )
returns the index of the net-flux for the fid flux.

This function searches for the flux fid in the original network, uses this to
resolve the WX-Pair, and returns the first index of the pair.

# NOTE
If the flux is NOT part of a WX-pair, then it returns the flux.
"""
function get_wxflux_net(net::Network , wt::WX_Transform , fid::String )
    idx_flux = find( x -> x==fid , net.fluxes )

    findpair = (wt.wparams .== idx_flux)
    swp      = sum(findpair[:])
    if(swp==1)
        pair = find( any(findpair,2) )[1]
        return wt.wparams[pair,1]
    else
        if(swp==0) # flux not in wx pair, return idx itself
            return idx_flux
        else
            error("Something wrong.. swp = $(swp)")
        end
    end
end

"""
  get_wxflux_xch(net::Network , wt::WX_Transform , fid::String )
returns the index of the net-flux for the fid flux.

This function searches for the flux fid in the original network, uses this to
resolve the WX-Pair, and returns the second index of the pair.

# NOTE
If the flux is NOT part of a WX-pair, then it returns -1
"""
function get_wxflux_xch(net::Network , wt::WX_Transform , fid::String )
    idx_flux = find( x -> x==fid , net.fluxes )

    findpair = (wt.wparams .== idx_flux)
    swp      = sum(findpair[:])
    if(swp==1)
        pair = find( any(findpair,2) )[1]
        return wt.wparams[pair,2]
    else
        if(swp==0) # flux not in wx pair, return idx itself
            return -1
        else
            error("Something wrong.. swp = $(swp)")
        end
    end
end

"""
  s_wx_flux(wxt::WX_Transform, flux::String)
returns whether given flux has wx parametrization or not
"""
function is_wx_flux(net::Network , wt::WX_Transform , fid::String )
    idx_flux = find( x -> x==fid , net.fluxes )
    if( any( wt.wparams[:] .== idx_flux ) )
        return true
    else
        return false
    end
end

"""
  find_wx_flux(wxt::WX_Transform, flux::String)
returns positions (flux[n],flux[x]) .
Note that either of the two original fluxes can be supplied here.
"""
function get_wx_flux_pair(net::Network , wt::WX_Transform , fid::String )
    idx_flux = find( x -> x==fid , net.fluxes )
    if( any( wt.wparams[:] .== idx_flux ) )
        find_pair = find( any( (wt.wparams .== idx_flux) , 2) )
        pair      = wt.wparams[find_pair,:]
        # sanity check
        if(pair[2]<pair[1]); error("Ordering of wx fluxes wrong.."); end
        return (pair[1],pair[2])
    else
        return (-1,-1)
    end
end

"""
  find_wx_flux(wxt::WX_Transform, flux::String)
returns the index of the given wx paramtrization parameter.
Notation is: flux_id[n] , flux_id[x] for net / xchange flux
"""
function find_wx_flux( net::Network, wxt::WX_Transform, flux::String)
    if(flux[end]==']' && flux[end-2]=='[')
        if( flux[end-1]=='n' || flux[end-1]=='N' )
            return get_wxflux_net(net::Network , wxt::WX_Transform , flux[1:end-3] )
        elseif( flux[end-1]=='x' || flux[end-1]=='X' )
            return get_wxflux_xch(net::Network , wxt::WX_Transform , flux[1:end-3] )
        end
    end
    # if no [x] / [n] is specified, return the flux index..
    return find(fi -> fi==flux ,net.fluxes)
end


"""
  create_all_wx_fluxnames(net::Network,wxt::WX_Transform)
creates all wx flux names
Notation is: flux_id[n] , flux_id[x] for net / xchange flux
"""
function create_all_wx_fluxnames(net::Network,wxt::WX_Transform)
    flux_indeces = collect(1:length(net.fluxes))

    test_wx_flux = (x::Int64) -> any( wxt.wparams[:] .== x )

    flux_names   = Vector{String}(length(net.fluxes))

    for zi in 1:length(net.fluxes)
        if(test_wx_flux(zi))
            wxp = get_wx_flux_pair(net,wxt,net.fluxes[zi])
            flux_names[ wxp[1] ] = string( net.fluxes[wxp[1]] , "[n]")
            flux_names[ wxp[2] ] = string( net.fluxes[wxp[1]] , "[x]")
        else
            flux_names[ zi ]     = net.fluxes[zi]
        end
    end
    return flux_names
end

"""
lc2polyhedron(net::Network , wx::WX_Transform , lcs::Vector{LinConstr})
returns Polyhedron
creates polyhedron from the given constraints

# NOTE: Format of  constraints: flux_id[x] is xchange flux, flux_id[n] is net flux
"""
function lc2polyhedron(net::Network , wxt::WX_Transform , lcs::Vector{LinConstr} )
    n = length(net.fluxes)

    # Init the flux resolution dictionary:
    # NOTE: for fluxes which are not wx-parametrized, and which therefore basically do not have
    # the [x]/[n] postfix, we map the [n] version to the non-postfix flux
    d_s::Dict{String,Int64}   = Dict{String,Int64}()
    #for zi=1:length(net.fluxes); d_s[net.fluxes[zi]] = zi; end
    wx_fluxes = create_all_wx_fluxnames(net,wxt)
    for zi in 1:length(wx_fluxes); print(string(zi," --> " , wx_fluxes[zi] ,"\n")) ;d_s[ wx_fluxes[zi] ] = zi; end
    # add the [n] aliases for non wx-parametrized fluxes
    for fi in wx_fluxes
        if(fi[end]!=']'); d_s[string(fi,"[n]")] = d_s[fi]; end
    end

    d_s["1"] = length(wx_fluxes)+1 # to collect the "1" data

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

"""
  get_N(net::Network,wxt::WX_Transform)
returns stoichiometric matrix for network with wx transform
"""
function get_N(net::Network,wxt::WX_Transform)
    # create N, then sort out [xch] fluxes..
    N = get_N(net)

    # remove all constraints involing "second" fluxes of wxparams
    for iwx in wxt.wparams[:,2]
        #icr = find( N[:,iwx] .!= 0 )
        #N = N[setdiff(collect(1:size(N,1)),icr),:]
        N[:,iwx] = 0 # this should be ok, because in every constraint of N, the forward flux should always be non-zero, if this flux is non-zero..
    end
    N_all_zero = find( all( N .==0 ,2) )
    N = N[ setdiff(collect(1:size(N,1)),N_all_zero),:]
    return N
end

"""
  get_flux_polyhedron(net::Network , wxt::WX_Transform , constraints::Vector{LinConstr} ; max_val=Inf )
creates the flux polyhedron for network with wx_transform and wx constraints
CONVENTION:
Here, all (net) fluxes are reversible by default
We limit the xch fluxes to the interval [0,1] by default.
"""
function get_flux_polyhedron(net::Network , wxt::WX_Transform , constraints::Vector{LinConstr} ; max_val=Inf )

     m_n = get_N(net,wxt)
     n   = size(m_n,2)
     p_0 = lc2polyhedron(net , wxt , constraints)

     d_s::Dict{String,Int64}   = Dict{String,Int64}()
     for zi=1:length(net.fluxes); d_s[net.fluxes[zi]] = zi; end

     # create [0,1] bounding box for xch fluxes
     ub =  Inf*ones(n)
     lb = -Inf*ones(n)
     lb[wxt.wparams[:,2]] = 0.0
     ub[wxt.wparams[:,2]] = 1.0

     # if max_val, create bouding box. [-max_val,max_val] around everything..
     if(max_val<Inf)
         ub[ub.>max_val] = max_val
         lb[ get_idx_sx_fluxes(net) ] = -max_val
     end

     p_1 =  Polyhedron( p_0 , bounding_box(lb,ub) )
     p_2 =  Polyhedron( p_1 , Polyhedron( zeros(0,n),zeros(0,1),m_n,zeros(size(m_n,1),1) ) )
     return p_2
end



struct optwx_parametrization
    rfd_p::rfd_parametrization
    wx_t::WX_Transform
end

function eval_opt2raw( tf::optwx_parametrization , x::Array{Float64,2})
    x_pwx   = tf.rfd_p.at_opt2raw(x)
    x_flux  = MetabolicFluxAnalysis.eval_wx_transform(tf.wx_t,x_pwx)
end

function get_optimal_wx_parametrization(net::Network,lcs::Vector{LinConstr},max_flux::Float64)
    (p_pwx,wtx)  = MetabolicFluxAnalysis.parametrization_wiechert_exchange( net, lcs, max_flux )

    # round wx_parameterspace polytope:
    rfd_parametrization_pwx = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( p_pwx ; remove_red_constr=1)

    # sample dat polytope unif:
    parametrization = optwx_parametrization(rfd_parametrization_pwx,wtx)
    return parametrization
end


# find flux  pairs with "exchange flux" ( <-> ), and reparametrize them by
# "net-flux" + "exchange flux".
# n: stoichiometric matrix
# sx_param optiosn:
# "log":  is "log4" , means, sxp(0) = 0.0001*max, sxp(1) = max
function parametrization_wiechert_exchange( net::MetabolicFluxAnalysis.Network , lcs::Vector{MetabolicFluxAnalysis.LinConstr} , max_flux::Float64 ; sx_param="log")#::Polyhedron )

    if(max_flux>1e6); print("Set max flux to 1e6.. need finite bound\n"); end
    max_flux = min(max_flux,1e6)

    p_constr = MetabolicFluxAnalysis.lc2polyhedron(net , lcs::Vector{LinConstr} )

    print("Wiechert SX Parametrization:\n")

    # find pairs in stoich. matrix:
    N = MetabolicFluxAnalysis.get_N(net)
    w_pairs = zeros(Int64,0,2)

    # store parametrized fluxes, to avoid problems with multiple reversible fluxes with same stoichiometry
    # e.g. to handle --> fum_a,fum_a_bw,fum_b,fum_b_bw
    parametrized_fluxes = Vector{Int64}()

    for zi in 1:size(N,2)
        for zj in zi+1:size(N,2)
            if( in(zi,parametrized_fluxes) || in(zj,parametrized_fluxes) ) # skip!
                continue;
            end

            if( rank(N[:,[zi,zj]])==1 )
                # and check if is in opposite direction:
                if( sign.(N[:,zi] ) == -sign.(N[:,zj]) )
                    print("pair: ($(zi),$(zj))  -   $(net.fluxes[zi]) , $(net.fluxes[zj]) ")

                    #check if there are other constraints on it:
                    cfree = !any( p_constr.A[:,[zi,zj]][:] .!=0 ) && !any( p_constr.G[:,[zi,zj]][:] .!=0 )
                    if(!cfree)
                        print(" --> skip parametrization, found in constraints..")
                    else
                        w_pairs = [w_pairs;[zi zj]]
                        push!( parametrized_fluxes , zi); push!( parametrized_fluxes , zj)
                        print(" --> parametrized!")
                    end
                    print("\n");
                end
                #w_pairs = [ w_pairs ; [zi zi]]
            end
        end
    end
    if( sx_param=="log" || sx_param=="log4")
        fxch = (xpxch::Float64) ->  exp( xpxch * (log(10000.))  ) * (max_flux/10000.)
    end
    if( sx_param=="log5")
        fxch = (xpxch::Float64) ->  exp( xpxch * (log(100000.))  ) * (max_flux/100000.)
    end
    if( sx_param=="log3")
        fxch = (xpxch::Float64) ->  exp( xpxch * (log(1000.))  ) * (max_flux/1000.)
    end
    if(sx_param=="abs")
        fxch = (xpxch::Float64) -> fxch * max_flux
    end

    ft = (xnet::Float64,xxch::Float64) -> (  (xnet>0) ? xnet+xxch : xxch , (xnet>0) ? xxch : -xnet + xxch )
    ftp = (xnet::Float64,pxxch::Float64) -> ft( xnet , fxch(pxxch) )

    wtx =WX_Transform(w_pairs,ftp)

    # create parameter space polytope:
    # add [0,1] bounds for all wx-parametrized sx fluxes
    # add [-max , max] for wx-parametrized net flux
    # and do the default thing: bound all (non-wx-parametrized) non-sx fluxes by [0,max_flux]
    # and all sx fluxes by [-max_flux,max_flux]
    # !! AND !!
    # Remove the sx parametrized flux from the stoichiometric matrix..
    p_v = get_flux_polyhedron_vanilla(net::Network , lcs::Vector{LinConstr} )

    # remove fluxes w_pairs[:,2] from the stoichiometric matrix:
    p_v.A[:,w_pairs[:,2]] = 0 #zeros(size(wx,1))


    # identify non-system-exchange fluxes
    fnsx = get_idx_non_sx_fluxes(net)
    #print(fnsx)

    lb =  Inf * -ones(length(net.fluxes))
    ub =  Inf *  ones(length(net.fluxes))
    lb[fnsx] = 0
    # if max_val, create bouding box. [0,max_val] for nsx fluxes, AND [-max_val,max_val] for sx fluxes
    if(max_flux<Inf)
        ub[ub.>max_flux] = max_flux
        lb[ get_idx_sx_fluxes(net) ] = -max_flux
    end
    lb[w_pairs[:,2]] = 0
    ub[w_pairs[:,2]] = 1
    lb[w_pairs[:,1]] = -max_flux
    ub[w_pairs[:,1]] = max_flux


    pp = Polyhedron( p_v , bounding_box(lb,ub) )

    return (pp,wtx)
end
