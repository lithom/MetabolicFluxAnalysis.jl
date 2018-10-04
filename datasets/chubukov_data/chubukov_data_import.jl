using ExcelFiles, DataFrames
using MetabolicFluxAnalysis

# conditions = excel sheet names
conditions = ["glucose","malate","malate-glucose","glycerol","pyruvate","gluconate","aspartate","glutamate succinate","fructose"]


map_fragments = Dict{String,Tuple{String,Vector{Integer}}}()

map_fragments["ALA-15"]  = ("Ala",[])
map_fragments["ALA-57"]  = ("Ala",[1,2,3])
map_fragments["ALA-85"]  = ("Ala",[1,3])

map_fragments["ASP-15"]  = ("Asp",[1,2,3,4])
map_fragments["ASP-57"]  = ("Asp",[1,2,3,4])
map_fragments["ASP-85"]  = ("Asp",[])
map_fragments["ASP159"] = ("Asp",[1,2,4])
map_fragments["ASP302"] = ("Asp",[1,2])

map_fragments["GLU-15"]  = ("Glu",[])
map_fragments["GLU-57"]  = ("Glu",[1,2,3,4,5])
map_fragments["GLU-85"]  = ("Glu",[1,2,4,5])
map_fragments["GLU159"]  = ("Glu",[1,2,4,5])
map_fragments["GLU302"]  = ("Glu",[1,2])

map_fragments["GLY-15"]  = ("Gly",[])
map_fragments["GLY-57"]  = ("Gly",[1,2])
map_fragments["GLY-85"]  = ("Gly",[2])

map_fragments["HIS-57"]  = ("His",[])
map_fragments["HIS159"] = ("His",[])

map_fragments["ILE-15"] = ("Ile",[])
map_fragments["ILE-85"] = ("Ile",[1,2,4,5,6])
map_fragments["ILE159"] = ("Ile",[1,2,4,5,6])

map_fragments["LEU-15"] = ("Leu",[])
map_fragments["LEU-85"] = ("Leu",[1,2,4,5,6])
map_fragments["LEU159"] = ("Leu",[1,2,4,5,6])

map_fragments["PHE-57"]  = ("Phe",[1,2,3,4,5,6,7,8,9])
map_fragments["PHE-85"]  = ("Phe",[1,2,3,4,5,6,7,9])
map_fragments["PHE159"]  = ("Phe",[1,2,3,4,5,6,7,9])
map_fragments["PHE302"]  = ("Phe",[])

map_fragments["PRO-85"]  = ("Pro",[])
map_fragments["PRO159"]  = ("Pro",[1,3,4,5])

map_fragments["SER-15"]  = ("Ser",[])
map_fragments["SER-57"]  = ("Ser",[1,2,3])
map_fragments["SER-85"]  = ("Ser",[1,3])
map_fragments["SER159"]  = ("Ser",[1,3])
map_fragments["SER302"]  = ("Ser",[])

map_fragments["THR-15"]  = ("Thr",[])
map_fragments["THR-57"]  = ("Thr",[1,2,3,4])
map_fragments["THR-85"]  = ("Thr",[])

map_fragments["TYR-15"]  = ("Tyr",[])
map_fragments["TYR-57"]  = ("Tyr",[])
map_fragments["TYR-85"]  = ("Tyr",[])
map_fragments["TYR159"]  = ("Tyr",[])

map_fragments["VAL-15"]  = ("Val",[])
map_fragments["VAL-57"]  = ("Val",[1,2,3,4,5])
map_fragments["VAL-85"]  = ("Val",[1,2,3,4])
map_fragments["VAL159"]  = ("Val",[1,2,3,4])
map_fragments["VAL302"]  = ("Val",[])






exp_m_glucose = Dict{Int64,String}()
exp_m_glucose[1] = "e1m1"
exp_m_glucose[2] = "e1m2"
exp_m_glucose[3] = "e2m1"
exp_m_glucose[4] = "e2m2"
exp_m_glucose[5] = "e2m3"
exp_m_mal = Dict{Int64,String}()
exp_m_mal[1] = "e1m1"
exp_m_mal[2] = "e1m2"
exp_m_mal[3] = "e2m1"
exp_m_mal[4] = "e2m2"
exp_m_mal[5] = "e2m3"
exp_m_malgluc = Dict{Int64,String}()
exp_m_malgluc[1] = "e1m1"
exp_m_malgluc[2] = "e1m2"
exp_m_malgluc[3] = "e2m1"
exp_m_malgluc[4] = "e2m2"
exp_m_malgluc[5] = "e2m3"
exp_m_glycerol = Dict{Int64,String}()
exp_m_glycerol[1] = "e1m1"
exp_m_glycerol[2] = "e1m2"
exp_m_glycerol[3] = "e2m1"
exp_m_glycerol[4] = "e2m2"
exp_m_glycerol[5] = "e2m3"
exp_m_pyruvate = Dict{Int64,String}()
exp_m_pyruvate[1] = "e1m1"
exp_m_pyruvate[2] = "e1m2"
exp_m_pyruvate[3] = "e2m1"
exp_m_pyruvate[4] = "e2m2"
exp_m_pyruvate[5] = "e2m3"
exp_m_pyruvate[6] = "e3m1"
exp_m_pyruvate[7] = "e3m2"
exp_m_pyruvate[8] = "e3m3"
exp_m_pyruvate[9] = "e4m1"
exp_m_pyruvate[10] = "e4m2"
exp_m_pyruvate[11] = "e4m3"
exp_m_gluconate = Dict{Int64,String}()
exp_m_gluconate[1] = "e1m1"
exp_m_gluconate[2] = "e1m2"
exp_m_gluconate[3] = "e2m1"
exp_m_gluconate[4] = "e2m2"
exp_m_gluconate[5] = "e2m3"
exp_m_aspartate = Dict{Int64,String}()
exp_m_aspartate[1] = "e1m1"
exp_m_aspartate[2] = "e1m2"
exp_m_aspartate[3] = "e1m3"
exp_m_aspartate[4] = "e2m1"
exp_m_aspartate[5] = "e2m2"
exp_m_aspartate[6] = "e2m3"
exp_m_glutamsuccin = Dict{Int64,String}()
exp_m_glutamsuccin[1] = "e1m1"
exp_m_glutamsuccin[2] = "e1m2"
exp_m_glutamsuccin[3] = "e1m3"
exp_m_glutamsuccin[4] = "e2m1"
exp_m_glutamsuccin[5] = "e2m2"
exp_m_glutamsuccin[6] = "e2m3"
exp_m_fructose = Dict{Int64,String}()
exp_m_fructose[1] = "e1m1"
exp_m_fructose[2] = "e1m2"
exp_m_fructose[3] = "e2m1"
exp_m_fructose[4] = "e2m2"
exp_m_fructose[5] = "e2m3"


conditions = ["glucose","malate","malate-glucose","glycerol","pyruvate","gluconate","aspartate","glutamate succinate","fructose"]

map_measurements = Dict{String,Dict{Int64,String}}()
map_measurements["glucose"] = exp_m_glucose
map_measurements["malate"] = exp_m_mal
map_measurements["malate-glucose"] = exp_m_malgluc
map_measurements["glycerol"] = exp_m_glycerol
map_measurements["pyruvate"] = exp_m_pyruvate
map_measurements["gluconate"] = exp_m_gluconate
map_measurements["aspartate"] = exp_m_aspartate
map_measurements["glutamate succinate"] = exp_m_glutamsuccin
map_measurements["fructose"] = exp_m_fructose



# parse everything:

measurements = Dict{ String, Vector{Tuple{SpeciesEMU,Tuple{Int64,Array{Float64}}}} }()

for zi in 1:length(conditions)

    condition_i = conditions[zi]
    measurements_i = map_measurements[ condition_i ]



    df = DataFrame(load("C:\\Users\\Thomas\\Desktop\\chubukov_data\\chubukov MDVs.xls", condition_i))

    # process the fragment mass description:
    list_of_fragments = Vector{Tuple{SpeciesEMU,Tuple{Int64,Array{Float64}}}}()

    for zj in 6:256
        str_frag = df[1][zj]
        str_frag_split = split(str_frag," ")
        frag = map_fragments[str_frag_split[1]]
        if(isempty(frag[2]))
            continue # fragment not clear..
        end

        str_mass = replace(str_frag_split[2],"m","")
        i_mass = parse(Int64,str_mass)

        m_values = Vector{Float64}()
        for ki in keys(measurements_i)
            # pos is (ki+1)
            push!( m_values , df[ki+1][zj] )
        end
        push!( list_of_fragments , ( constructSpeciesEMU(frag[1],frag[2]) , (i_mass,m_values) ) )
    end

    measurements[condition_i] = list_of_fragments
end

v1 = Vector([1,2,3])
