using MetabolicFluxAnalysis

# implements the volume computation algorithm presented in
# I. Z. Emiris and V. Fisikopoulos , 2014 , Efﬁcient random-walk methods for approximating polytope volume.
# currently, it uses slower (uniform direction) hr sampling instead of the
# proposed chr sampling. This will be implemented soon.
function vol_mmc_reversed( G , h , num_samples , n_steps_initial , n_steps_rep , rho )

    d = size(G,2)

    # num iterations to do, such that rho^m >= d^d
    m = Integer( ceil( ( d * log(d) ) / log(rho) ) )

    # construct the radii that we have to do in the rounded polytope..
    # i.e. all from 1.0 to
    r_m = rho.^( (0:m) / d )

    print("r_m: $(r_m)")

    # round polytope:
    if(true)
        #tic()
        # compute MVE
        # E = {L^(-1)*x + f: norm(x) <= 1}
        (f,L) = mve_run(G, h)
        # rescale to isotropic ball
        h = h - G*f
        G = (G/L)/sqrt(d+2)
        #time_rounding = toc()
    end

    # Generate initial samples:
    #tic()
    # warm-up samples..
    print("\nWarmup.. \n")
    #switch(hr_conf.method)
    #    case {'HitAndRun'}
    X0 = sample_unif_hr( MetabolicFluxAnalysis.Polyhedron(G,h[:,:]) , zeros(d,num_samples) ,  n_steps_initial )
    #    case 'CoordHR',
            #%X0 = fullCoordHR(G,h,100*r_m(end),zeros(d,1),num_samples,n_steps_rep);
            #X0 = fullCoordHR(G,h,100*r_m(end),zeros(d,1),num_samples,n_steps_initial_decorr);
    #end
    #time_initialsamples = toc()


    ## Run over the reverse mmc iterations..
    X = X0

    factors = ones(Integer(m)) #Vector{Float64}(m)
    n_out_stored = zeros(Integer(m))#Vector{Int64}(m)
    n_resampled_tot = 0

    #tic()
    # initial radius:
    for zi=m:-1:2
        # ball radius of iteration:
        r_old = sqrt(d+2) * r_m[zi]
        r_new = sqrt(d+2) * r_m[zi-1]

        idx_out = sqrt.( sum( X.^2 , 1 ) ) .> r_new
        n_out   = sum(idx_out)
        ratio   = num_samples ./ (num_samples-sum(idx_out))

        n_out_stored[zi] = n_out
        # store result
        factors[zi] = (num_samples/(num_samples-n_out))
        n_resampled_tot = n_resampled_tot + n_out
        # Remove samples that are out:
        X(:,idx_out) = []

        # Resample:
        if( (n_out>0) && (zi > 2) )
            # select random samples which will be decorrelated..
            idx_rep = randsample( size(X,2) , n_out ) # should work here
            X_rep = X[ : , idx_rep ]

            # Run Hit-and-run steps:
            #switch(hr_conf.method),
            #    case {'HitAndRun','CoordinateHR'},
            #X_rep = XUniformHitAndRun_SingleCore( X_rep , struct('A',G,'b',h) , struct('L',eye(d),'f',zeros(d,1),'r',r_new,'L_is_eye',1) , n_steps_rep , hr_conf , hr_output_conf  );
            #X_rep = hrstep!( G , h , r_new , X_rep , n_steps_rep)
            X_rep = hrstep_cchr!!( G , h , r_new , X_rep , n_steps_rep)
            #    case 'MexCoordHR',
            #        X_rep = fullCoordHR(G,h,r_new,X_rep,1,n_steps_rep);
            #end
            #X[:,end+1:end+n_out] = X_rep;
            X = [ X  X_rep ]
        end
        # it/allit   factora factorb   tot_resampled
        #fprintf(' %04d / %04d   %d/%d    %d    \n',   m-zi+1,m  , num_samples,num_samples-n_out , n_resampled_tot );

        #if(total_debug), plotcolors = linspecer(m); plot(X(1,1:500),X(2,1:500),'.','Color',plotcolors(zi,:)); end
    end
    #time_mmc = toc()

    # original volume of ellipsoid corresopnding to the smallest ball:
    #lvol = (d/2)*log(pi) - gammaln(d/2+1) - sum(log(diag(L)))
    lvol = (d/2)*log(pi) - lgamma(d/2+1) - sum(log.(diag(L)))
    #%lvol = (d/2)*log(pi) - gammaln(d/2); %- sum(log(diag(L)));

    # multiply with the estimated factors..
    lvol = lvol + sum( log.(factors[2: length(factors) ] ) )

    logvol = lvol
    vol    = exp(lvol)
    #MoreInfo.time_Rounding       = time_rounding;
    #MoreInfo.time_InitialSamples = time_initialsamples;
    #MoreInfo.time_MMC            = time_mmc;
    #MoreInfo.time_Total          = toc(tic_start);

    #MoreInfo.factors             = factors(2:end);
    #MoreInfo.n_out               = n_out_stored(2:end);
    #MoreInfo.logvol_mve          = (d/2)*log(pi) - gammaln(d/2+1) - sum(log(diag(L)));
end
