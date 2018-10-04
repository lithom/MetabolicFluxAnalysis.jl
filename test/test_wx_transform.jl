

using MetabolicFluxAnalysis
using CoordinateTransformations

#reload("MetabolicFluxAnalysis")

# current tests..
fluxi_1_f = open("testdata/test_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)
fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )


## Try wiechert sx parametrization:
(p_pwx,wxt)  = MetabolicFluxAnalysis.parametrization_wiechert_exchange( fluxi1.net, fluxi1.linear_constraints["ConstraintsExchange"] , 2000. )

MetabolicFluxAnalysis.create_all_wx_fluxnames(fluxi1.net,wxt)

# round wx_parameterspace polytope:
rfd_parametrization_pwx = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( p_pwx ; remove_red_constr=1)

# sample dat polytope unif:
x_opt_unif = MetabolicFluxAnalysis.sample_unif_hr(rfd_parametrization_pwx.p_opt,rfd_parametrization_pwx.cheby_p_opt,1000,1000)


x_pwx_unif = rfd_parametrization_pwx.at_opt2raw(x_opt_unif)
x_flux     = MetabolicFluxAnalysis.eval_wx_transform(wxt,x_pwx_unif)
#x_flux     = @enter MetabolicFluxAnalysis.eval_wx_transform(wtx,x_pwx_unif)

x_flux[:,1:10]

# check that it is steady state fluxes:
@test all( ( ( MetabolicFluxAnalysis.get_N(fluxi1.net) * x_flux ).^2 )[:] .< 1e-7 )



## test with bigger network:

# load file
#fluxi_1_f = open("datasets/complete_mfa/cmfa_01.fluxi")
fluxi_1_f = open("testdata/cmfa_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)
fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )

constraints_config_01 = ["Constraints_FixedRatios","Constraints_BiomassRatios","Constraints_Measurements_01","Constraints_IrreversibilityPublication"]
constraints_01 = reduce( (x,y)->[x;y] , map( x -> fluxi1.linear_constraints[x] , constraints_config_01 ) )

fp1 = MetabolicFluxAnalysis.get_flux_polyhedron( fluxi1.net , constraints_01 , max_val=2000 )
wxp = MetabolicFluxAnalysis.get_optimal_wx_parametrization(fluxi1.net,constraints_01,2000.)


# sample:
x_wx_opt =  MetabolicFluxAnalysis.sample_unif_hr( wxp.rfd_p.p_opt , wxp.rfd_p.cheby_p_opt , 5000 , 200 )
x_flux    =  MetabolicFluxAnalysis.eval_opt2raw( wxp , x_wx_opt )
