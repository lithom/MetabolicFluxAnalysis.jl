
#@testset "Fluxi.Basics_01" begin

    using MetabolicFluxAnalysis
    lca = MetabolicFluxAnalysis.parseLC(" 3.2 A + 1   B + C + 2   D  ")


    s1 = MetabolicFluxAnalysis.Species("s1",4)
    s2 = MetabolicFluxAnalysis.Species("s2",0)
    s3 = MetabolicFluxAnalysis.Species("a2",3)
    ss = [s1;s2;s3]
    #print( sort(ss) )


    sas1 = MetabolicFluxAnalysis.SAS("sa",1)
    sas2 = MetabolicFluxAnalysis.SAS("sa",2)
    sas3 = MetabolicFluxAnalysis.SAS("sb",1)
    sas4 = MetabolicFluxAnalysis.SAS("sa",2)

    isequal(sas2,sas4)

    hash(sas2)==hash(sas4)
    ss1 = Set( [sas1;sas2;sas3;sas4] )
    @test length(ss1)==3


    sas1a = MetabolicFluxAnalysis.SAS("sa",1)
    sas2a = MetabolicFluxAnalysis.SAS("sa",2)
    sas1b = MetabolicFluxAnalysis.SAS("sa",1)
    sas2b = MetabolicFluxAnalysis.SAS("sa",2)

    ss1 = Set( [sas1a;sas2a] )
    ss2 = Set( [sas1b;sas2b] )
    @test ss1==ss2



    ras1 = MetabolicFluxAnalysis.RAS("f1",MetabolicFluxAnalysis.SAS("sa",1),1)
    ras2 = MetabolicFluxAnalysis.RAS("f1",MetabolicFluxAnalysis.SAS("sb",1),1)
    ras3 = MetabolicFluxAnalysis.RAS("f1",MetabolicFluxAnalysis.SAS("sa",1),1)
    ps1 = Set([ras1;ras2;ras3])
    @test length(ps1)==2


    pas1 = MetabolicFluxAnalysis.PAS("f1",MetabolicFluxAnalysis.SAS("sa",1),1)
    pas2 = MetabolicFluxAnalysis.PAS("f1",MetabolicFluxAnalysis.SAS("sb",1),1)
    pas3 = MetabolicFluxAnalysis.PAS("f1",MetabolicFluxAnalysis.SAS("sa",1),1)
    ps1 = Set([pas1;pas2;pas3])
    @test length(ps1)==2




    m1 = MetabolicFluxAnalysis.ASTEMUMassRatio( MetabolicFluxAnalysis.SpeciesEMU( [MetabolicFluxAnalysis.SAS("A",1);MetabolicFluxAnalysis.SAS("A",3)] ) , 1.0 )
    m2 = MetabolicFluxAnalysis.ASTEMUMassRatio( MetabolicFluxAnalysis.SpeciesEMU( [MetabolicFluxAnalysis.SAS("A",1);MetabolicFluxAnalysis.SAS("A",3)] ) , 1.0 )
    @test hash(m1) == hash(m2)
    @test isequal( m1.semu , m2.semu )
    @test isequal(m1,m2)

#end
