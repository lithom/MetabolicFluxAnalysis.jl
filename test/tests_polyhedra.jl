
#reload("MetabolicFluxAnalysis.jl")
using MetabolicFluxAnalysis
using CoordinateTransformations

# current tests..
fluxi_1_f = open("testdata/test_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)
fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )


# now get and then parametrize the flux polytope
fp = MetabolicFluxAnalysis.get_flux_polyhedron( fluxi1.net , fluxi1.linear_constraints["ConstraintsExchange"] , max_val=200 )

# compute transform to make full-dimensional:
( fp_fd , at_fd ) = MetabolicFluxAnalysis.find_minimal_affine_subspace(fp)
fp_fd_a = MetabolicFluxAnalysis.remove_redundant_constraints( fp_fd )
( fp_rounded , at_rounded ) = MetabolicFluxAnalysis.round_polytope(fp_fd_a)

at_opt2flux = CoordinateTransformations.compose(at_fd , at_rounded)



# test chebyshev comp.
x_int_cheby = MetabolicFluxAnalysis.solve_chebyshev(fp_fd.G,fp_fd.h[:])

@test all( fp_fd.G * x_int_cheby .< fp_fd.h )
print( fp_fd.h - fp_fd.G * x_int_cheby  )

## test with toy polytopes..
d_test = 100
g_test = [eye(d_test);-eye(d_test)]; h_test = 10*rand(2*d_test,1);
a_test = rand(2,d_test)       ; b_test = 0.1 * rand(2,1);
f_test = MetabolicFluxAnalysis.Polyhedron(g_test,h_test,a_test,b_test)
rfd_parametrization = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( f_test ; remove_red_constr=1)
@assert MetabolicFluxAnalysis.isinside( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt[:,:] )

@test MetabolicFluxAnalysis.isinside( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt[:,:] )

# generate some unifrom samples over polytope:
xx = MetabolicFluxAnalysis.sample_unif_hr( rfd_parametrization.p_opt , rfd_parametrization.cheby_p_opt , 1000 ,1000 );

xx_raw = MetabolicFluxAnalysis.eval_opt2raw(rfd_parametrization,xx)

@test MetabolicFluxAnalysis.checksamples( rfd_parametrization.p_opt , xx     ; eq_tol=1e-4)
@test MetabolicFluxAnalysis.checksamples( rfd_parametrization.p_raw , xx_raw ; eq_tol=1e-4)





##


# some more basic tests
# check that things behave correctly for unbounded polyhedra
using MetabolicFluxAnalysis
bb_01 = MetabolicFluxAnalysis.bounding_box([-10.;-10;-10],[10.;10;20])
unbounded_polyhedron = MetabolicFluxAnalysis.Polyhedron( bb_01.G[1:1,:],bb_01.h[1:1,:],bb_01.A[:,:],bb_01.b[:,:] )
lpm_result = MetabolicFluxAnalysis.solve_lp_max_v(unbounded_polyhedron,[0.;1;0])




#(mve_f,mve_E) =  MetabolicFluxAnalysis.mve_run(fp_fd_a.G,fp_fd_a.h[:])





# test rounding transform:
cube_g = [1.0*eye(4);-eye(4)]
cube_h = 1.0*collect(1:8)

p_c = MetabolicFluxAnalysis.Polyhedron(cube_g,cube_h[:,:])

# remove_redundant_constraintss
param1= fp_fd_a.G
param2 = fp_fd_a.h[:]
param3 = [79.514 ;52.79]
(mve_f,mve_E) =  MetabolicFluxAnalysis.mve_solver( param1,param2,param3 )


(mve_f,mve_E) =  MetabolicFluxAnalysis.mve_run(cube_g,cube_h)


# transform somehow:

p_fixed = MetabolicFluxAnalysis.Polyhedron(fp_fd,MetabolicFluxAnalysis.bounding_box([-1000.;-1000],[1000.;1000]) )

(mve_f,mve_E) = MetabolicFluxAnalysis.mve_run( p_fixed.G , p_fixed.h[:] )
(mve_f,mve_E) = MetabolicFluxAnalysis.mve_run( fp_fd.G , fp_fd.h[:] )

# compute rounding transform:
(mve_f,mve_E) = MetabolicFluxAnalysis.mve_run( fp_fd.G , fp_fd.h[:] )
