using MetabolicFluxAnalysis


x0  = 0.1*rand(2,1000)
tic()
MetabolicFluxAnalysis.hrstep!( [eye(2);-eye(2)] , [ones(2);ones(2)] , 0.8*sqrt(2) , x0 , 2000)
toc()


x0  = 0.1*rand(2,1000)
tic()
MetabolicFluxAnalysis.hrstep_cchr!( [eye(2);-eye(2)] , [ones(2);ones(2)] , 0.8*sqrt(2) , x0 , 2000)
toc()


function benchmark_pencilbox( d::Int , n_samples::Int , n_steps::Int )
    A = 1.0 * [eye(d);-eye(d)]
    b = 1.0 * [1:d;zeros(d)]

    x0  = 0.1*rand(d,n_samples)
    r   = 10.0

    tic()
    MetabolicFluxAnalysis.hrstep_cchr!( A , b , r , x0 , n_steps)
    t1 = toc()
    tic()
    MetabolicFluxAnalysis.hrstep!(      A , b , r , x0 , n_steps)
    t2 = toc()

    return (t1 , t2)
end


benchmark_pencilbox( 100,1000,1000 )


#using Plots
#gr()
#plot()
#scatter!(x0[1,:],x0[2,:])
#plot!()
(lambda_min,lambda_max,x_intersection_min,x_intersection_max) = MetabolicFluxAnalysis.xrayshoot_centerball(1.2,0.0*rand(2,200),rand(2,200))
#using Plots
#gr()
#plot()
#scatter!(x_intersection_max[1,:],x_intersection_max[2,:])
#plot!()



using MetabolicFluxAnalysis
pa1 = [eye(2);-eye(2)]
pb1 = [1.;2;3;4] # size  should be (1+3)*(2+4) = 24
MetabolicFluxAnalysis.vol_mmc_reversed(   pa1 , pb1 , 500,200,10,2  )
MetabolicFluxAnalysis.vol_mmc_reversed(   pa1 * qr(rand(2,2))[1] , pb1 , 500,200,10,2  )
MetabolicFluxAnalysis.vol_mmc_reversed(   pa1 * qr(rand(2,2))[1] , pb1 , 500,200,10,2  )
