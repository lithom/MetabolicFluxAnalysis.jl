import Base.isless
import Base.isequal
import Base.==
import Base.hash

import Base.show

"""
SAS has fields species::String and pos::Int
it represents an atom slot in a species
"""
struct SAS
    species::String
    pos::Int64
end

function Base.isless(x::SAS,y::SAS)
    if( x.species < y.species ); return true; end
    if( x.species == y.species)
        return x.pos < y.pos
    else
        return false
    end
end
function ==(x::SAS,y::SAS)
    return isequal(x.species,y.species) && x.pos==y.pos
end
function Base.hash(x::SAS)
    return hash(x.species)+hash("xxx$(x.pos)")
end
function Base.isequal(x::SAS,y::SAS)
    return isequal(x.species,y.species) && x.pos==y.pos
end


"""
RAS has fields flux::String, sas::SAS and rpos::Int
it represents an atom slot in a reactand molecule of a specific flux
"""
struct RAS
    flux::String
    sas::SAS
    rpos::Int64
end
# sort: 1. flux 2. species, 3. (!!) rpos 4. sas.pos
function Base.isless(x::RAS,y::RAS)
    if( x.flux < y.flux ); return true; end
    if( x.flux == y.flux )
        if( x.sas.species < y.sas.species ); return true; end
        if( x.species == y.species)
            if(x.rpos < y.rpos)
                return true
            end
            if(x.rpos == y.rpos)
                return x.sas.pos < y.sas.pos
            end
        end
        return false
    end
    return false
end
function ==(x::RAS,y::RAS)
    return  x.flux==y.flux && isequal(x.sas,y.sas) && x.rpos==y.rpos
end
function Base.hash(x::RAS)
    return hash("$(x.flux)ab$(x.sas.species)xx$(x.rpos)yy$(x.sas.pos)")
end
function Base.isequal(x::RAS,y::RAS)
    return x.flux==y.flux && Base.isequal(x.sas,y.sas) && x.rpos==y.rpos
end
function Base.show(io::IO, x::RAS)
  print(io, "RAS:($(x.flux),$(x.sas.species),$(x.sas.pos));rpos=$(x.rpos))")
end

"""
PAS has fields flux::String, sas::SAS and rpos::Int
it represents an atom slot in a product molecule of a specific flux
"""
struct PAS
    flux::String
    sas::SAS
    ppos::Int64
end
# sort: 1. flux, 2. species, 3. (!!) rpos 4. sas.pos
function Base.isless(x::PAS,y::PAS)
    if(x.flux <  y.flux); return true; end
    if(x.flux == y.flux)
        if( x.sas.species < y.sas.species ); return true; end
        if( x.sas.species == y.sas.species)
            if(x.ppos < y.ppos)
                return true
            end
            if(x.ppos == y.ppos)
                return x.sas.pos < y.sas.pos
            end
        end
        return false
    end
    return false
end

function ==(x::PAS,y::PAS)
    return x.flux==y.flux && Base.isequal(x.sas,y.sas) && x.ppos==y.ppos
end
function Base.hash(x::PAS)
    return return hash("$(x.flux)ab$(x.sas.species)xx$(x.ppos)yy$(x.sas.pos)")
end
function Base.isequal(x::PAS,y::PAS)
    return x.flux==y.flux && Base.isequal(x.sas,y.sas) && x.ppos==y.ppos
end
function Base.show(io::IO, x::PAS)
  print(io, "PAS:($(x.flux),$(x.sas.species),$(x.sas.pos));ppos=$(x.ppos))")
end

"""
MassDistributionMeasurement has fields sass::Set{SAS} and data::Dict{Int64,Array{Float64,1}}
It represents a mass distribution measurement for the set of SAS given by .sass
.data contains measurement data for the different isotopic weights, i.e. the
key is the weight, and the value is the data.
Format of .data[m] is usually either:
  .data[m][1] : mean of measurement
  .data[m][2] : std. dev. of measurement
or (in case of multiple data points):
  .data[m][i] : i'th measurement
"""
struct MassDistributionMeasurement
    sass::Set{SAS}
    data::Dict{Int64,Array{Float64,1}}
end


function MassDistributionMeasurement( sass::Set{SAS} , data::Vector{Tuple{Int64,Array{Float64}}} )
    dict_data = Dict{Int64,Array{Float64,1}}()
    for zi in 1:length(data)
        dict_data[data[zi][1]] = data[zi][2]
    end
    for zm in 0:length(sass)
        if( ! contains(==,collect(keys(dict_data)),zm) )
            dict_data[zm] = [0.0;0.0;0.0]
        end
    end
    return MassDistributionMeasurement(sass,dict_data)
end

#struct ConfiguredSAS
#    sass::Set{SAS}
"""
ConfiguredSEMU contains fields .semu::Set{SAS} and .config::Vector{Int64}
.config specifies unlabeled / labeled atoms (0/1)
"""
struct ConfiguredSEMU
    semu::Set{SAS}
    config::Vector{Int64}
end
function ==(x::ConfiguredSEMU,y::ConfiguredSEMU)
    return Base.isequal(x.semu,y.semu) && x.config==y.config
end
function Base.hash(x::ConfiguredSEMU)
    return hash(x.semu)+hash(x.config)
end
function Base.isequal(x::ConfiguredSEMU,y::ConfiguredSEMU)
    return Base.isequal(x.semu,y.semu) && x.config==y.config
end
function Base.show(io::IO, x::SAS)
    print(io, "SAS:($(x.species),$(x.pos))")
end

ConfiguredSEMUMixture = Dict{ConfiguredSEMU,Float64}

SpeciesEMU = Set{SAS}
ReactandEMU = Set{RAS}
ProductEMU = Set{PAS}


function constructSpeciesEMU(species::String,pos::Vector{Integer})
    return constructSpeciesEMU(species,convert(Array{Int64},pos))
end

function constructSpeciesEMU(species::String,pos::Vector{Int64})
    semu1 = SpeciesEMU()
    for zi=1:length(pos)
        push!(semu1,SAS(species,pos[zi]))
    end
    return semu1
end


function Base.isless(a::SpeciesEMU,b::SpeciesEMU)
    if(length(a)<length(b)); return true; end
    if(length(a)==length(b))
        as = sort(collect(a))
        bs = sort(collect(b))
        for i in 1:length(a)
            if(as[i]<bs[i])
                return true
            else
                if(bs[i]<as[i])
                    return false
                end
            end
        end
        return false
    end
    return false
end

function Base.show(io::IO, x::SpeciesEMU)
    xs = sort(collect(x))
    if( all( map( x1 -> x1.species == xs[1].species , xs ) ) )
        print(io,"SEMU::$(xs[1].species)::(")
        str_slots = ""
        for xi in xs; str_slots = string(str_slots,"$(xi.pos),") end
        print(io,str_slots[1:length(str_slots)-1])
        print(io,")")
    else
        print(io,"SEMU::")
        for xi in xs
            print(io,"($(xi.species),$(xi.pos))")
        end
        #print(io,"")
    end
end

function get_positions_of_semu(s::SpeciesEMU)
    if(!isinsinglespecies(s)); error("getpositionofsemu requires single species EMU!"); end
    return map( si->si.pos  ,s)
end

function remu2semu(re::ReactandEMU)
    return SpeciesEMU( map( x -> x.sas , re ) )
end
function pemu2semu(pe::ProductEMU)
    return SpeciesEMU( map( x -> x.sas , pe ) )
end


function isinsinglespecies(se::SpeciesEMU)
    if(isempty(se)); return true; end
    e_species = collect(se)[1].species
    return all( map( x -> x.species==e_species , collect(se) ) )
end




"""""
matching_mass_combinations( target_mass::Int64 , max_individual_masses::Vector{Int64} )
return Vector{Vector{Int64}}

returns all combinations of masses which (i) sum up to target_mass, and where
mass i has value inbetween 0 and max_individual_masses[i]

Example:
matching_mass_combinations( 2 , [1,3,2] )
returns: [0,0,2] , [0,1,1] , [0,2,0] , [1,0,1] , [1,1,0]

matching_mass_combinations( 3 , [1,4] )
return [0,4] , [1,3]
etc..
"""
function matching_mass_combinations( target_mass::Int64 , max_individual_masses::Vector{Int64} )
    # compute recursively: take the first mim value, and loop over all possibilities..
    combis = Vector{Vector{Int64}}()

    if(length(max_individual_masses)==1)
        # end of recursion..
        if(max_individual_masses[1]>=target_mass)
            return [[target_mass]]
        else
            return Array{Array{Int64,1},1}(0)
        end
    end

    mim_rec = max_individual_masses[2:end]
    for mi_w in 0:max_individual_masses[1]
        if(target_mass-mi_w<0); continue; end
        mmc_i = matching_mass_combinations( target_mass-mi_w , mim_rec )
        mmc_i_2 = map( xi -> [mi_w;xi] , mmc_i )
        append!( combis , mmc_i_2 )
    end
    return combis
end




function benchmarks_01()

    f_randi = (max_i) -> trunc(Int,rand()*max_i*0.99999)+1

    flux_names    = map( x -> string("v_",randstring()) ,1:100)
    species_names = map( x -> string("s_",randstring()) ,1:100)

    sas1::Vector{SAS} = Vector{SAS}()
    for zi=1:200
        push!( sas1 , SAS(species_names[f_randi(100)] , f_randi(12) ) )
    end

    semus1::Vector{SpeciesEMU} = Vector{SpeciesEMU}()
    for zi=1:200
        i_spc   = species_names[f_randi(100)]

        i_slots_pre = collect( 1:f_randi(12) )
        i_slots     = shuffle(i_slots_pre)
        i_sl        = i_slots[1:f_randi(length(i_slots))]
        sas_i::Vector{SAS} = Vector{SAS}()
        for zj=1:length(i_sl); push!(sas_i,SAS(i_spc,i_sl[zj])); end
        semui       = SpeciesEMU(Set(sas_i))
        push!( semus1 , semui )
    end

    # test dicts:
    tic()
    d_semu::Dict{SpeciesEMU , Float64} = Dict{SpeciesEMU , Float64}()
    for zi=1:1e4
        d_semu = Dict{SpeciesEMU , Float64}()
        for zj=1:100
            d_semu[semus1[zj]] = rand()
        end
    end
    t_dicts_1 = toc()

    # test dicts 2:
    semus_to_choose = collect( keys( d_semu ) )

    tic()
    sum_1::Float64 = 0.0
    for zi=1:1000000
        sum_1 += 0.0001 * d_semu[ semus_to_choose[ f_randi(length(semus_to_choose)) ] ]
    end
    t_dicts_2 = toc()
    println("\ndict acccess speed: $(1000000/t_dicts_2) / s  \n")
end
