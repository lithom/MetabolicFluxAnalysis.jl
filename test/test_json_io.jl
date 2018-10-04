
##
reload("MetabolicFluxAnalysis")
using MetabolicFluxAnalysis


# load file
#fluxi_1_f = open("datasets/antoniewicz_toy_network/toy_network_01.fluxi")
fluxi_1_f = open("testdata/toy_network_01.fluxi")
fluxi_1_s = readstring(fluxi_1_f)
fluxi1 = MetabolicFluxAnalysis.parseFluxi( fluxi_1_s )


#
str_json_01 = MetabolicFluxAnalysis.export_json(fluxi1)

print(str_json_01)

# try import..
fluxi_imported = MetabolicFluxAnalysis.import_json(str_json_01)


fluxi_1_reconstructed = MetabolicFluxAnalysis.import_json_to_fluxi(str_json_01)
