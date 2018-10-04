using JuMP
using Clp

using ASTInterpreter2

#
# All credits for this piece of software to Yin Zhang
# 

function solve_chebyshev(A::Array{Float64,2},b::Array{Float64,1})
# Use a hyperplane representation of the polyhedron
# i.e. P = { x | Ax ≤ b }
# m = 4
#A = [ 2. 1 ; 2 -1 ; -1 2 ; -1 -2]
#b = ones(4)
(sm,sn) = size(A)

# Build JuMP model
m = Model(solver=ClpSolver())
@variable(m, r)
@variable(m, x[1:sn])
@objective(m, Max, r)
sA = sparse(A)
@constraint( m , sA * x + sqrt.(sum(sA.^2,2)) * r .<= b )
#@addConstraint(m, P[i=1:4], dot(A[i],x) + norm(A[i])*r <= b[i])

# Solve with any available LP solver
solve(m)

# Where is the center?
#print(getvalue(x))
return getvalue(x)

end

#%  MVE_RUN  Find the maximum volume ellipsoid of a polytope.
#%
#%     [f,E] = MVE_RUN(A, b) finds the maximum volume inscribed ellipsoid of
#%     the polytope {x: A*x <= b}. The polytope is assumed to be full-
#%     volume. The ellipse is returned as the pair (f,E) such that
#%         Ell = {E^(-1)*x + f: norm(x) <= 1}
#%     or equivalently as
#%         Ell = {x:  norm(E*(x-f)) <= 1},
#%     The matrix E is upper triangular.
#%
#%     NOTE: An ellipsoid can be represented as:
#%     (1) {x: (x-f)'*(E'*E)*(x-f) <= 1}  , or
#%     (2) {x: {E^(-1)*x + f, norm(x) <= 1}
#%
#%
#%     [f,E] = MVE_RUN(A, b, x0) finds the maximum volume inscribed
#%     ellipsoid using the starting point x0. It must be an interior point,
#%     i.e. b - A*x0 > 0.
#
#%--------------------------------------
#% Yin Zhang, Rice University, 07/29/02
#%--------------------------------------
#% Markus Uhr, ETH Zurich, 10/04/11
#%--------------------------------------
 # Thomas Liphardt, June 2018
 #--------------------------------------

function mve_run(A::Array{Float64,2},b::Array{Float64,1} ; x_internal=[] )
    tol = 1.e-8
    maxiter = 200

    x0 = NaN * ones(size(A,2))
    if(isempty(x_internal))
        x0 = solve_chebyshev(A,b)
    else
        x0 = x_internal
    end
    (f,E,msg,y,z,iter) = mve_solver( A , b , x0 ; maxiter=maxiter , tol=tol )
    return(f,E)
end

#function [x,E,msg,y,z,iter] = mve_solver(A,b,x0,maxiter,tol)
#  MVE_SOLVER  Driver for solving the maximum volume ellipsoid problem
#%
#%  See also: MVE_RUN
#
#%    {v:  v = x + E^(-1)*s, ||s|| <= 1}
#%  inscribing a full-dimensional polytope
#%          {v:  Av <= b}
#%  Input:  A, b --- defining the polytope
#%          x0 --- interior point (Ax0 < b)
#%  Output: x --- center of the ellipsoid
#%          E2 --- E'*E
#%--------------------------------------
#% Yin Zhang, Rice University, 07/29/02
#%--------------------------------------
#% Markus Uhr, ETH Zurich, 10/04/11
#%--------------------------------------
#% Thomas Liphardt, June 2018
#%--------------------------------------


function mve_solver(A_p::Array{Float64,2},b::Array{Float64,1},x0::Array{Float64,1}; maxiter::Int64=50 , tol::Float64=1.e-4 )
    (m, n) = size(A_p);
    bnrm = norm(b);

    #if ~exist('maxiter', 'var')
    #    maxiter = 50;
    #end;
    #if ~exist('tol', 'var')
    #    tol = 1.e-4;
    #end;

    minmu::Float64 = 1.e-8
    tau0::Float64 = .75

    bmAx::Array{Float64,1} = b - A_p*x0
    if any(bmAx .<= 0)
        error("x0 not interior")
    end

    A = sparse(1:m,1:m,1./bmAx)*A_p
    b = ones(m)

    x = zeros(n)
    y = ones(m)
    bmAx = b

    z::Array{Float64,1} = NaN * ones(n)

    R::Array{Float64,2} = NaN * ones(n,n)

    print("\n  Residuals:   Primal     Dual    Duality  logdet(E)\n");
    print("  --------------------------------------------------\n");

    msg = 0

    Adx::Array{Float64,1} = NaN * ones(n)# init Adx
    astep::Float64 = NaN # init astep

    iter::Int64 = 0

    #%----- loop starts -----
    for iter = 1:maxiter

        if (iter > 1)
            bmAx = bmAx .- (astep*Adx) # TODO check if this is the same..
        end

        Y = sparse(1:m,1:m,y)
        #E2 = inv(full(transpose(A)*Y*A));
        #Q = A*E2*A';
        furz = 7
        #R = chol( transpose(A) *Y*A)
        R = chol( Hermitian( transpose(A) *Y*A ) )
        Q = A*(R\( transpose(R) \ transpose(A) ))


        h = sqrt.(diag(Q))
        if (iter == 1)
            t = minimum(bmAx./h)
            y = y/t^2
            h = t*h
            z = bmAx-h
            z[z.<1.e-1] = 1.e-1
            Q = t^2*Q
            Y = Y/t^2
        end
        yz = y.*z
        yh = y.*h
        gap = sum(yz)/m
        rmu = minimum([.5; gap])*gap
        rmu = maximum([rmu;minmu])

        R1 = - transpose(A) * yh
        R2 = bmAx - h - z
        R3 = rmu - yz

        r1 = norm(R1,Inf)
        r2 = norm(R2,Inf)
        r3 = norm(R3,Inf)
        res = maximum([r1;r2;r3])
    #%     objval = log(det(E2))/2;
        #objval = full(-sum(log.(diag(R))))
        objval = (-sum(log.(diag(R))))

        @printf("  iter %3i  ", iter);
        #fprintf('%9.1e %9.1e %9.1e  %9.3e\n', r2,r1,r3,objval);
        @printf("%9.1e %9.1e %9.1e  %9.3e\n", r2,r1,r3,objval);
        if ( res < tol*(1+bnrm) && rmu <= minmu )
            #fprintf('  Converged!\n');
            x = x + x0
            msg=1
            break
        end

        #print("\nCHECK 01\n");

        YQ = Y*Q;
        YQQY = YQ.*YQ'
        y2h = 2*yh
        YA = Y*A
        G  = YQQY + sparse(1:m,1:m,max.(1.e-12,y2h.*z))
        T = G \ (sparse(1:m,1:m,h+z)*YA)
        ATP = (sparse(1:m,1:m,y2h)*T-YA)'

        #print("\nCHECK 02\n");

        R3Dy = R3./y;
        R23 = R2 - R3Dy;
        dx = (ATP*A)\(R1 + ATP*R23)
        Adx = A*dx

        dyDy = G\(y2h.*(Adx - R23))
        dy = y.*dyDy
        dz = R3Dy - z.*dyDy

        ax = -1/minimum([-Adx./bmAx; -.5])
        ay = -1/minimum([ dyDy; -.5])
        az = -1/minimum([dz./z; -.5])
        tau = maximum( [tau0; 1 - res])
        astep = tau*minimum([1;ax;ay;az])

        #print("\nCHECK 03\n");

        x = x + astep*dx
        y = y + astep*dy
        z = z + astep*dz
    end

    E = R

    return (x,E,msg,y,z,iter)
end

##
