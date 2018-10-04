using Clp
using JuMP
using CoordinateTransformations

import Base.copy


# redefine CoordinateTransformations map evaluation funtion (is this allowed? is this bad?)
function (trans::CoordinateTransformations.AffineMap{M, V})(x::Array{Float64,2}) where {M, V}
    broadcast( + , trans.m * x , trans.v )
end


mutable struct Polyhedron
    G::Array{Float64,2}
    h::Array{Float64,2}
    A::Array{Float64,2}
    b::Array{Float64,2}
end

function Polyhedron(G::Array{Float64,2},h::Array{Float64,2})
    return Polyhedron(G,h,Array{Float64,2}(0,size(G,2)),Array{Float64,2}(0,1))
end


# intersection of two polyhedra
function Polyhedron(pa::Polyhedron,pb::Polyhedron)
    return Polyhedron( [pa.G;pb.G] , [pa.h;pb.h] , [pa.A;pb.A] , [pa.b;pb.b] )
end

function Base.copy(p::Polyhedron)
    return Polyhedron(copy(p.G),copy(p.h),copy(p.A),copy(p.b))
end


function isinside(p::Polyhedron , x::Array{Float64,2} ; eq_tol=1e-4)
    dx =  broadcast( - , p.h , p.G * x )

    if( ~isempty(p.A) )
        eq_viol = broadcast( - , p.A * x , b )
        return all( dx[:] .> 0 ) && all(eq_viol .< eq_tol)
    else
        return all( dx[:] .> 0 )
    end
end

function checksamples(p::Polyhedron , x::Array{Float64,2} ; eq_tol=1e-4)
    dx =  broadcast( - , p.h , p.G * x )
    if( ~isempty(p.A) )
        eq_viol = broadcast( - , p.A * x , p.b )
        print("Max Eq. Viol:   $(maximum(eq_viol[:]))")
        print("Min Ineq. Dist: $(minimum(dx[:]))")
        return all( dx[:] .> 0 ) && all(eq_viol .< eq_tol)
    else
        print("Min Ineq. Dist: $(minimum(dx[:]))")
        return all( dx[:] .> 0 )
    end

end



function bounding_box_old(lbs::Vector{Float64},ubs::Vector{Float64})
    n = length(lbs)
    idx_lb_not_inf = find( x -> ( !isinf(x) )  , lbs )
    idx_ub_not_inf = find( x -> ( !isinf(x) )  , ubs )
    G::Array{Float64,2} = Array{Float64}(0,length(lbs))
    h::Array{Float64,2} = Array{Float64}(0,1)
    for zi=1:length(idx_lb_not_inf)
        zii = idx_lb_not_inf[zi]
        gi = zeros(1,n); gi[zii] = -1;
        G = [G ;  gi ]
        h = [h ; -lbs[zii] ]
    end
    for zi=1:length(idx_ub_not_inf)
        zii = idx_ub_not_inf[zi]
        gi = zeros(1,n); gi[zii] = 1;
        G = [G ;  gi ]
        h = [h ; ubs[zii] ]
    end
    return Polyhedron( G,h )
end


function bounding_box(lbs::Vector{Float64},ubs::Vector{Float64})
    n = length(lbs)
    idx_lb_not_inf = find( x -> ( !isinf(x) )  , lbs )
    idx_ub_not_inf = find( x -> ( !isinf(x) )  , ubs )
    G_s = spzeros(0,length(lbs))
    h::Array{Float64,2} = Array{Float64}(0,1)
    for zi=1:length(idx_lb_not_inf)
        zii = idx_lb_not_inf[zi]
        gi = spzeros(1,n); gi[zii] = -1;
        G_s = [G_s ;  gi ]
        h = [h ; -lbs[zii] ]
    end
    for zi=1:length(idx_ub_not_inf)
        zii = idx_ub_not_inf[zi]
        gi = spzeros(1,n); gi[zii] = 1;
        G_s = [G_s ;  gi ]
        h = [h ; ubs[zii] ]
    end
    return Polyhedron( full(G_s),h )
end




#function find_minimal_affine_subspace(p::Polyhedron ; min_diam=1e-6 )
function find_minimal_affine_subspace( p::Polyhedron )

    conf_min_diam = 1e-6
    output_level = 2
    nc = size(p.G,1)
    d  = size(p.G,2)

    # sort out "zero constraints.."
    idx_zero_constr = []
    for zi in 1:nc
        if(norm(p.G[zi,:]) < 1e-6 )
            push!(idx_zero_constr,zi)
        end
    end
    if(!isempty(idx_zero_constr))
        print(" WARNING! $(length(idx_zero_constr)) infinite inequality constraint found.. ")
        p.G = p.G[ setdiff( (1:nc),(idx_zero_constr)) , : ]
        p.h = p.h[ setdiff( (1:nc),(idx_zero_constr)) , : ]
        nc = size(p.G,1)
    end

    ## Part 1: Check all inequality constraints
    fmins = Vector{Float64}(nc)
    fmaxs = Vector{Float64}(nc)

    x_boundary_min = Array{Float64,2}(d,nc)
    x_boundary_max = Array{Float64,2}(d,nc)

    #print( size(p.G) , size(p.h) , size(p.A) , size(p.b) )
    if(output_level>=1)
        print("Probe $(d) inequality constraints: ")
    end
    for zi in 1:nc
        if(output_level==1); if( mod( ceil(nc/8) )==0 ); print("."); end; end
        if(output_level==2); print("\n$(@sprintf("%3d",zi)) : "); end

        #(lp_a_f,lp_b_f,lp_a_x,lp_b_x) = solve_lp_minmax(p,p.G[zi,:])
        #print("  " ,norm(p.G[zi,:]), "  ")
        (lp_a_f,lp_b_f,lp_a_x,lp_b_x) = solve_lp_minmax(p, p.G[zi,:] / norm(p.G[zi,:]) )
        fmins[zi] = lp_a_f; fmaxs[zi,1] = lp_b_f
        x_boundary_min[:,zi] = lp_a_x
        x_boundary_max[:,zi] = lp_b_x

        if(output_level==2); print("diam= $(@sprintf("%10.6f",lp_b_f-lp_a_f))"); end
    end
    if(output_level>=1); print("\n"); end

    x_boundary = [x_boundary_min  x_boundary_max]
    f_diams = (fmaxs-fmins)
    c_ineq = f_diams .>= conf_min_diam
    c_ineq_to_eq = .~ c_ineq


    if(output_level>2)
        #println("\nDistances: ")
        #println(f_diams)
    end

    if(output_level>0)
        println( "\nMin considered Diameter: $(minimum(f_diams[c_ineq] ))\n" )
        if(any(c_ineq_to_eq))
            println(   "Max non-considered Diameter: $(maximum(f_diams[c_ineq_to_eq] ))\n" )
        else
            println(   "All dimensions considered, polyhedron is full-dimensional")
        end
    end

    ## Part 2: Create the new polytope
    G = p.G[c_ineq,:]; h = p.h[c_ineq,1]
    A = [ p.A ; p.G[c_ineq_to_eq,:] ]
    b = [ p.b ; p.h[c_ineq_to_eq,:] ]

    ## PART 3: Compute nullspace and reparametrize:
    x0_ip     = sum(x_boundary,2) / (2*nc) # just take the mean of all computed points..
    null_A    = nullspace(A)

    # reparametrize:
    G_1 = G * null_A
    h_1 = h - ( G * x0_ip )
    # That should be it.. (?)
    p_new = Polyhedron(G_1,h_1)

#    return p_new

    at::AffineMap = AffineMap(null_A,x0_ip)
    return p_new , at
    # x0 = zeros(size(G_1,2),1) # origin is in reparametrized polytope.
end


function round_polytope(p::Polyhedron)
    if(!isempty(p.A)); error("Rounding works only for polytopes in ineq. representation!"); end

    (f,E) = mve_run(p.G,p.h[:])

    G_r = p.G * inv(E)
    h_r = p.h - p.G * f
    p_r = Polyhedron(G_r,h_r)

    at::AffineMap = AffineMap(inv(E),f)
    return p_r , at
end

# p = Polyhedron( ones(1,4) , ones(4,1) , [eye(4);-eye(4)] , 2*ones(4,1) )

function solve_lp_minmax(p::Polyhedron,v::Vector{Float64})
    m = Model(solver=ClpSolver())
    @variable(m, x[1:length(v)])
    #@objective(m, Max, AffExpr( x , v' , 0 ) )
    @objective(m, Min, v'*x )
    #@constraint( m , sparse(p.A) * x .<= p.b )
    sG = sparse(p.G)
    sA = sparse(p.A)
    @constraint( m , sG * x .<= p.h )
    @constraint( m , sA * x .== p.b )
    # Solve with any available LP solver
    solve(m)
    lp_a_x = getvalue(x)
    lp_a_f = getobjectivevalue(m)

    @objective( m , Max , v'*x )
    solve(m)
    lp_b_x = getvalue(x)
    lp_b_f = getobjectivevalue(m)


    return (lp_a_f,lp_b_f,lp_a_x,lp_b_x)
end


function solve_lp_max_v(p::Polyhedron,v::Vector{Float64})
    m = Model(solver=ClpSolver())
    @variable(m, x[1:length(v)])
    #@objective(m, Max, AffExpr( x , v' , 0 ) )
    @objective(m, Max, v'*x )
    @constraint( m , p.A * x .<= p.b )
    @constraint( m , p.G * x .<= p.h )
    @constraint( m , p.A * x .== p.b )
    # Solve with any available LP solver
    solve(m)
    lp_a_x = getvalue(x)
    lp_a_f = getobjectivevalue(m)
    return (lp_a_f,lp_a_x)
end

function remove_redundant_constraints(p::Polyhedron)
    red = map( xi -> test_redundancy(p,xi) , 1:size(p.G,1) )
    p2 = Polyhedron( copy( p.G[ .~ red , : ] ) , copy(p.h[ .~ red , : ])  , copy(p.A) , copy(p.b) )
    return p2
end

function test_redundancy(p::Polyhedron,i::Int64)
    p2 = copy(p)
    p_g = p2.G[i,:]
    p_h = p2.h[i,1]
    p2.G = p2.G[ (1:size(p.G,1)).!=i,:]
    p2.h = p2.h[ (1:size(p.G,1)).!=i,:]

    (x_f , x_x) = solve_lp_max_v(p2,p_g)

    if( isnan(x_f) )
        # whatever happend here.. play it safe and say we need it.
        # currently, the solve_lp_max_v function returns NaN for unbounded, so this is important..
        return false
    end

    # SAFETY CUSHION!!
    # to be on the safe side wrt. removing redundant constraints..
    CONF_redundance_test_cushion = 0.001;

    if( x_f > p_h-CONF_redundance_test_cushion)
        return false
    else
        return true
    end
 end


"""
heuristic way of identifying the constraints which make a polyhedron empty.
works by relaxing all constraints, until polyhedron is non-empty
"""
function heuristically_find_infeasible_constraints( p::Polyhedron ; relaxation_step = 0.1 , max_tries=500 )

    G_A_p = [p.A;-p.A]
    h_b_p = [p.b;-p.b]

    P_original = Polyhedron( [p.G;G_A_p] , [p.h;h_b_p] )
    Pi         = copy(P_original)

    test_polytope = (p::Polyhedron) -> !isnan( solve_lp_max_v(p,ones(size(p.G,2)))[1] )

    cnt = 0
    print("Relax:\n")
    polytope_ok = test_polytope(Pi)
    while(cnt<max_tries && ! polytope_ok)
        cnt += 1
        print("."); if(mod(cnt,60)==0);print("\n");end;
        Pi.h[:] = Pi.h + relaxation_step
        polytope_ok = test_polytope(Pi)
    end
    print("\nSucces! Polytope feasible.. Now shrink back:\n")
    # now check for every constraint, if the shrunk back version is feasible / infeasible..
    idx_problems = Vector{Int64}()

    for zi=1:length(Pi.h)
        P2 = copy(Pi)
        P2.h[zi] = P_original.h[zi]
        if(!test_polytope(P2))
            push!(idx_problems,zi)
            print("x"); if(mod(zi,60)==0);print("\n");end;
        else
            print("."); if(mod(zi,60)==0);print("\n");end;
        end
    end

    # and remap indeces:
    problems_eq   = filter( x->x >  size(p.G,1) , idx_problems) .- size(p.G,1)
    problems_ineq = filter( x->x <= size(p.G,1) , idx_problems)

    return (problems_ineq,problems_eq)
end
