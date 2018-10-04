using MetabolicFluxAnalysis


# return values: [ lambda_min , lambda_max , x_intersection_min , x_intersection_max ]
# xrayshoot considers the system of linear inequalities A*x<b, and assumes
# that all columns of x0 lie inside the solution set.
# Then this function returns the most negative / most postitive lambdas
# such that x0_i + lambda_i * v_i <= A*b holds.
#
#
#TEST:
#A1=[eye(2);-eye(2)]; b1=[1;1;1;1];
# lambda_min , lambda_max , x_intersection_min , x_intersection_max = XRayShoot(A1,b1,zeros(2,4),[[1;1],[1;0],[-1;1],[0.1;0.2]]);
# lambda_min , lambda_max , x_intersection_min , x_intersection_max = XRayShoot(A1,b1,zeros(2,1),[[1;1],[1;0],[-1;1],[0.1;0.2]]);
function xrayshoot( A::Array{Float64,2} , b::Array{Float64,1} , x0::Array{Float64,2} , v::Array{Float64,2} )
    # if only one point, then replicate it once for every supplied
    # direction
    if ( (size(x0,2)==1) && (size(v,2)>1) )
        x0 = repmat(x0,1,size(v,2))
    end

    ## 1. compute intersections with constraints:
    AZ = A*v
    S::Array{Float64,2} = broadcast(-,b,A*x0)./AZ
    smax = copy(S)
    smax[AZ.<=0] = Inf
    smax = minimum(smax,1)
    lambda_max = smax;

    smin = S
    smin[AZ.>=0] = -Inf
    #print(smin)
    smin = maximum(smin,1)
    lambda_min = smin

    ## also return the intersecting points:
    #if(nargout >= 3),
    x_intersection_min = x0 + broadcast( * , lambda_min , v )
    x_intersection_max = x0 + broadcast( * , lambda_max , v )

    return lambda_min , lambda_max , x_intersection_min , x_intersection_max
end


# hr step for each col of x will be performed in replace
# random numbers are generated by the Julia GLOBAL_RNG
function hrstep!(A::Array{Float64,2} , b::Array{Float64,1} , x::Array{Float64,2} , nsteps::Int64)
    (d,n) = size(x)
    for zi=1:nsteps
        rv = randn(d,n)
        #rv = broadcast( / , rv , sum(rv.^2,1) )
        ru = rand(1, n)
        (lambda_min , lambda_max , x_intersection_min , x_intersection_max) = xrayshoot( A , b , x , rv )
        x[:,:] = x_intersection_min + broadcast( * , x_intersection_max-x_intersection_min , ru )
    end
    #return xi
end



function xrayshoot_centerball(r::Float64,x0::Array{Float64,2},v::Array{Float64,2})
    n=size(x0,2)
    #compute  for (x'=x+t*v)
    # (x^2)=r^2
    # t^2*v^2+t*(2xv)+x^2-r^2=0
    ea = sum( v.*v , 1); eb= 2*sum(v.*x0,1); ec= sum(x0.*x0,1)-r^2;
    eD = eb.^2 - 4*ea.*ec
    if( any( eD.<0 )); print("eD<0 -> kaputt\n"); end

    ex = [ (-eb+sqrt.(eD))./(2*ea) ; (-eb-sqrt.(eD))./(2*ea) ]
    ex_g0 = ex[1,:].>0
    lambda_max = NaN * ones(1,n)
    lambda_min = NaN * ones(1,n)
    lambda_max[1,ex_g0]  = ex[1,ex_g0]
    lambda_max[1,.~ex_g0] = ex[2,.~ex_g0]
    lambda_min[1,ex_g0]  = ex[2,ex_g0]
    lambda_min[1,.~ex_g0] = ex[1,.~ex_g0]

    x_intersection_min  = x0 + broadcast(*,v,lambda_min)
    x_intersection_max  = x0 + broadcast(*,v,lambda_max)
    return (lambda_min,lambda_max,x_intersection_min,x_intersection_max)
end


"""
  hrstep!(A::Array{Float64,2} , b::Array{Float64,1} , r::Float64 , x::Array{Float64,2} , nsteps::Int64)

hr step with uniform direction choice for each col of x will be performed inside
the intersection of the polytope G*x<h and the ball with center origin and
radius r.
NOTE: random numbers are generated by the Julia GLOBAL_RNG
"""
function hrstep!(A::Array{Float64,2} , b::Array{Float64,1} , r::Float64 , x::Array{Float64,2} , nsteps::Int64)
    #TODO implement..

    (d,n) = size(x)
    for zi=1:nsteps
        rv = randn(d,n)
        #rv = broadcast( / , rv , sum(rv.^2,1) )
        ru = rand(1, n)

        (lambda_min , lambda_max , x_intersection_min , x_intersection_max) = xrayshoot( A , b , x , rv )
        # and
        (lambda_min_cb , lambda_max_cb , x_intersection_min_cb , x_intersection_max_cb) = xrayshoot_centerball( r , x , rv )

        ix_min_cb = lambda_min_cb .> lambda_min
        x_intersection_min[ : , ix_min_cb[:] ] = x_intersection_min_cb[:,ix_min_cb[:]]
        ix_max_cb = lambda_max_cb .< lambda_max
        x_intersection_max[ : , ix_max_cb[:] ] = x_intersection_max_cb[:,ix_max_cb[:]]

        x[:,:] = x_intersection_min + broadcast( * , x_intersection_max-x_intersection_min , ru )
    end
    #return xi
end

"""
  hrstep_cchr!(A::Array{Float64,2} , b::Array{Float64,1} , r::Float64 , x::Array{Float64,2} , nsteps::Int64)

hr step with cchr direction choice for each col of x will be performed inside the
intersection of the polytope G*x<h and the ball with center origin and radius r
NOTE: random numbers are generated by the Julia GLOBAL_RNG
"""
function hrstep_cchr!(A::Array{Float64,2} , b::Array{Float64,1} , r::Float64 , x::Array{Float64,2} , nsteps::Int64)
    #TODO implement..

    d  = size(A,2)
    m  = size(A,1)

    ns = size(x,2)

    rand_dirs  = rand( 1:d , nsteps , ns)
    rand_steps = rand( nsteps , ns)

    for zi in 1:size(x,2)
        xi = x[:,zi]
        mx_b_minus_ax = b - A * xi

        for zs in 1:nsteps
            rdi::Int = rand_dirs[zs,zi]

            #sb::Float64 =  Inf
            #sa::Float64 = -Inf

            # test lin. constraints
            #ti = mx_b_minus_ax ./ A[:,rdi]
            #ti_idx_pos = ti .> 0.
            #s_a::Float64 = maximum( [ -Inf ; ti[.~ti_idx_pos] ] )
            #s_b::Float64 = minimum( [  Inf ; ti[  ti_idx_pos] ] )
            #s_a::Float64 = maximum( ti[.~ti_idx_pos] )
            #s_b::Float64 = minimum( ti[  ti_idx_pos] )
            s_a::Float64 = -Inf
            s_b::Float64 =  Inf
            for zm in 1:m
                #if( ti[zm] > 0.0 )
                ti_m = mx_b_minus_ax[zm] / A[zm,rdi]
                if( ti_m > 0.0 )
                    s_b = min(s_b,ti_m)
                else
                    s_a = max(s_a,ti_m)
                end
            end


            # test radius:
            quad_a::Float64 = 1.0
            quad_b::Float64 = 2*xi[rdi]
            quad_c::Float64 = -(r*r) + sum( xi.*xi )
            quad_D::Float64 = quad_b*quad_b - 4*quad_a*quad_c

            quad_q = -(quad_b + (quad_b>0?1.0:-1.0) * sqrt(quad_D))/ 2.0;
            b_xa::Float64 = quad_q/quad_a;
            b_xb::Float64 = quad_c/quad_q;

            if(b_xa<b_xb)
                s_a = (b_xa>s_a)?b_xa:s_a;
                s_b = (b_xb<s_b)?b_xb:s_b;

            else
                s_a = (b_xb>s_a)?b_xb:s_a;
                s_b = (b_xa<s_b)?b_xa:s_b;
            end

            step_bounds_gap::Float64 = 0.000001
            t_i = s_a +  (s_b-s_a) *  step_bounds_gap + (1.0-2*step_bounds_gap) * rand_steps[zs,zi]
            #t_i = s_a +  (s_b-s_a) * rand_steps[zs,zi]

            xi[rdi] = xi[rdi] + t_i
            mx_b_minus_ax -= A[:,rdi] * t_i
        end
        x[:,zi] = xi
    end
end

function sample_unif_hr( p::Polyhedron , x::Array{Float64,1} , nsamples::Int64 , nsteps::Int64 )
    x2 = repmat(x,1,nsamples)
    hrstep!(p.G , p.h[:] , x2 , nsteps)
    return x2
end
function sample_unif_hr( p::Polyhedron , x::Array{Float64,2} , nsteps::Int64 )
    x2 = x #repmat(x,1,nsamples)
    hrstep!(p.G , p.h[:] , x2 , nsteps)
    return x2
end