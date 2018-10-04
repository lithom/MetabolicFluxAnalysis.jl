
using MetabolicFluxAnalysis

s1 = MetabolicFluxAnalysis.Species("s1",4)
s2 = MetabolicFluxAnalysis.Species("s2",0)
s3 = MetabolicFluxAnalysis.Species("s3",3)

# test parseReactand:
# "sA()" --> (1.0,"sA","")
# " 4.0 sB " --> (4.0,"sB","")
# " sC(a-b-X) " --> (1.0,"sC","a-b-X")
 (stoich1,sid1,atm1) = MetabolicFluxAnalysis.parseReactand(" sA() ")
 (stoich2,sid2,atm2) = MetabolicFluxAnalysis.parseReactand(" 4.0 sB ")
 (stoich3,sid3,atm3) = MetabolicFluxAnalysis.parseReactand("  sC(a-b-X) ")

@test stoich1==1.0
@test sid1=="sA"
@test atm1==""
@test stoich2==4.0
@test sid2=="sB"
@test atm2==""
@test stoich3==1.0
@test sid3=="sC"
@test atm3=="a-b-X"




f1 = MetabolicFluxAnalysis.parseFlux( "TestFlux", " A(1-2-3) + B(A-b-x) + B(f-g-j) + X(u) --> C(1-2-3-x-b) + D(A) + D(f) + X(g-j-u) " )
#MetabolicFluxAnalysis.sas(f1)

## benchmark basic functions
#function create_rand_network(n_fluxes::Int64, n_species::Int64 , avg_spec_per_flux:: )
