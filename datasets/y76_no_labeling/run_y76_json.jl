
import JSON

using MetabolicFluxAnalysis

# load file
json_1_f = open("datasets/y76_no_labeling/y76Model.json")
json_1_s = readstring(json_1_f)

json_data = JSON.parse(json_1_s)

imp_s_result = MetabolicFluxAnalysis.import_json_species( json_data["metabolites"] )
imp_f_result = MetabolicFluxAnalysis.import_json_fluxes( json_data["reactions"] )

import_a = MetabolicFluxAnalysis.import_json( json_data )

net_a = import_a[1]
constraints_a = import_a[2]

tic()
fp  = MetabolicFluxAnalysis.get_flux_polyhedron(net_a , constraints_a ; max_val=Inf)
toc()
sx_fluxes = MetabolicFluxAnalysis.get_sx_fluxes(net_a)

fva_result = MetabolicFluxAnalysis.fva(net_a,constraints_a,collect(1:100))
fva_lb     = fva_result[1]
fva_ub     = fva_result[2]
#rfp = MetabolicFluxAnalysis.get_rounded_fulldimensional_parametrization( fp ; remove_red_constr=1)
