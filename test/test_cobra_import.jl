

#using Requests
#ecoliModel = Requests.get("http://gcrg.ucsd.edu/sites/default/files/Attachments/Images/downloads/Ecoli_core/ecoli_core_model.mat")
#save(ecoliModel, "ecoli_core_model.mat")

n_cobra = MetabolicFluxAnalysis.import_cobra("testdata/ecoli_core_model.mat","model")
#(myN ,myN_p , myN_n) = getN( n_cobra )
(myN ,myN_p , myN_n) = MetabolicFluxAnalysis.get_N( n_cobra ; sparse_n=true )

# benchmark stoich. matrix generation..:
sns = Array{Any,1}(10_000)
tic()
for zi=1:10_000
    ( sns_i , snxx ,  snyy ) = MetabolicFluxAnalysis.get_N(n_cobra, sparse_n=rand()>0.5 )
    sns[zi] = sns_i
end
toc()
