using MetabolicFluxAnalysis
# test parser:

fluxi_1_f = open("testdata/test_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)



MetabolicFluxAnalysis.parseFluxi(fluxi_1_s)
s1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )

@test length( s1.net.species ) == 6
@test length( s1.net.fluxes )  == 9
