using MetabolicFluxAnalysis

reload("MetabolicFluxAnalysis")

#s_data = open("D:\\Thomas\\Downloads\\IsoDesign-v1.2.1\\influx_s-v2.11.1\\test\\e_coli.ftbl")
s_data = open("testdata/e_coli.ftbl")


ftbl_data = MetabolicFluxAnalysis.parse_ftbl_data(s_data)


fluxi_from_ftbl = MetabolicFluxAnalysis.ftbl_to_fluxi(ftbl_data)

MetabolicFluxAnalysis.get_wxflux_net(fluxi_from_ftbl.fluxidata.net,fluxi_from_ftbl.wx_transform,"tk1")
MetabolicFluxAnalysis.get_wxflux_xch(fluxi_from_ftbl.fluxidata.net,fluxi_from_ftbl.wx_transform,"tk1")
MetabolicFluxAnalysis.find_wx_flux(fluxi_from_ftbl.fluxidata.net,fluxi_from_ftbl.wx_transform,"tk1")

MetabolicFluxAnalysis.create_all_wx_fluxnames(fluxi_from_ftbl.fluxidata.net,fluxi_from_ftbl.wx_transform)[16]

all_wx_fluxes = MetabolicFluxAnalysis.create_all_wx_fluxnames(fluxi_from_ftbl.fluxidata.net,fluxi_from_ftbl.wx_transform)

# create wx parametrization-based flux polytope..
wx_constraints = [ fluxi_from_ftbl.wx_constraints_fluxes ; fluxi_from_ftbl.wx_constraints_eq_and_ineq ; fluxi_from_ftbl.wx_constraints_flux_measurements ]

for wci in wx_constraints
    println(wci)
end



#fp_wx_constraints = MetabolicFluxAnalysis.lc2polyhedron(fluxi_from_ftbl.fluxidata.net,fluxi_from_ftbl.wx_transform, wx_constraints[1:10] )
#fp_wx_bounded = MetabolicFluxAnalysis.Polyhedron( fp_wx_constraints , MetabolicFluxAnalysis.bounding_box(zeros(170),2*ones(170)) )
#fp_wx_rp = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( fp_wx_bounded ; remove_red_constr=1)


fp_wx = MetabolicFluxAnalysis.get_flux_polyhedron( fluxi_from_ftbl.fluxidata.net , fluxi_from_ftbl.wx_transform , wx_constraints ; max_val=20. )


for zi=1:132
    vi = 1.0*(1:132 .==zi )  #rand(132)
    result_i = MetabolicFluxAnalysis.solve_lp_minmax(fp_wx,vi./norm(vi))
    print("$(zi) : $(all_wx_fluxes[zi]) : $(result_i[1]) - $(result_i[2]) \n")
end

# TODO: make this work reliably, i.e. figure out what the problems are during parametrization..
if(false)
    fp_wx_rp = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( MetabolicFluxAnalysis.Polyhedron(fp_wx,MetabolicFluxAnalysis.bounding_box( -100*ones(132) , 100*ones(132) )  ) ; remove_red_constr=0)

    vix = rand(132,5000)
    vix = broadcast( / , vix , sqrt.(sum(vix,1)) )
    results_lp = MetabolicFluxAnalysis.solve_lp_minmax.([fp_wx],  [ slicedim( vix,2,i) for i in 1:1000 ] )
    minimum( map( ri -> ri[2] , results_lp) )

    for zi=1:10000
        vi = rand(132)
        result_i = MetabolicFluxAnalysis.solve_lp_minmax(fp_wx,vi./norm(vi))
        print("$(vi[1]) - $(vi[2]) \n")
    end
end

## analyze fluxes..







##





MetabolicFluxAnalysis.heuristically_find_infeasible_constraints( fp_wx ; relaxation_step = 0.1 , max_tries=500 )

find( fp_wx.A[35,:].!=0)
find( fp_wx.A[134-size(fp_wx.A,1),:].!=0)
find( fp_wx.A[141-size(fp_wx.A,1),:].!=0)




fp_wx.G[1,14]
fp_wx.h[1]



find( fp_wx.G[108,:].!=0 )

find( x -> x>2 , sum( [fp_wx.G;fp_wx.A] .!=0 , 1 ) )


find( x -> x>2 , sum( [fp_wx.G] .!=0 , 1 ) )
