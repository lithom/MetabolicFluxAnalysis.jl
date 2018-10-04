
"""
FluxiData has fields:
  net::Network
  measurements::Dict{String,Vector{MassDistributionMeasurement}}
  substrate_configs::Dict{String,Vector{ConfiguredSEMUMixture}}
  linear_constraints::Dict{String,Vector{LinConstr
"""
struct FluxiData
    net::Network
    measurements::Dict{String,Vector{MassDistributionMeasurement}}
    substrate_configs::Dict{String,Vector{ConfiguredSEMUMixture}} # 
    linear_constraints::Dict{String,Vector{LinConstr}}
end


"""
parseFluxi(data::String)
returns FluxiData object

parses fluxi input
"""
function parseFluxi(data::String)

    ## Parse species:
    species::Vector{Species} = Vector{Species}()
    # 1: Find Species Section:
    d_species = extractSection( data , "Species" )

    # 2. split and loop over lines:
    s_lines = split(d_species,"\n")
    for ss in s_lines
        ss::String = preprocessLine(ss)
        if(isempty(ss)); continue; end
        ss_s = split(ss,":")
        push!(species,Species(strip(ss_s[1]),parse(Int64,strip(ss_s[2]))) )
    end

    ## Parse fluxes:
    fluxes::Vector{Flux} = Vector{Flux}()
    # 1: Find Species Section:
    d_fluxes = extractSection( data , "Fluxes" )

    # 2. split and loop over lines:
    f_lines = split(d_fluxes,"\n")
    for fs in f_lines
        fs::String = preprocessLine(fs)
        if(isempty(fs)); continue; end
        fs_s = split(fs,":")
        push!(fluxes, parseFlux( strip(String(fs_s[1])) , strip(String(fs_s[2])) ) )
    end

    ## create network from species and fluxes:
    fluxes_str::Vector{String}        = Vector{String}()
    fluxes_db::Dict{String,Flux}     = Dict{String,Flux}()
    species_str::Vector{String}       = Vector{String}()
    species_db::Dict{String,Species} = Dict{String,Species}()
    for fi in fluxes
        push!(fluxes_str,fi.id); fluxes_db[fi.id] = fi;
    end
    for si in species
        push!(species_str,si.id); species_db[si.id] = si;
    end
    net::Network = Network(fluxes_str,species_str,fluxes_db,species_db)


    ## Parse all measurements sections
    d_ms = extractNamedSections(String(data), "MeasurementsSet")
    dict_measurements = Dict{String,Vector{MassDistributionMeasurement}}()
    for msname in keys(d_ms)
        str_msi = d_ms[msname]
        mdml::Vector{MassDistributionMeasurement} = parseMeasurementSet(str_msi)
        dict_measurements[msname] = mdml
    end

    ## Parse all substrate config sections
    d_scs = extractNamedSections(String(data), "SubstrateConfiguration")
    #print(d_scs)
    dict_scs = Dict{String,Vector{ConfiguredSEMUMixture}}()
    for scname in keys(d_scs)
        str_scsi = d_scs[scname]
        scsi::Vector{ConfiguredSEMUMixture} = parseSubstrateMixture(str_scsi)
        dict_scs[scname] = scsi
    end

    ## Parse all linear constraints sections
    d_lcs = extractNamedSections(String(data), "LinearConstraintsSet")
    dict_lcs = Dict{String,Vector{LinConstr}}()
    for lcname in keys(d_lcs)
        str_lcsi = d_lcs[lcname]
        lcsi::Vector{LinConstr} = parseLinConstr(str_lcsi)
        dict_lcs[lcname] = lcsi
    end

    return FluxiData(net,dict_measurements,dict_scs,dict_lcs)
end


function extractSection(data::String, name::String)
    r_find_s = "$(name)Start:"
    r_find_e = "$(name)End:"
    mm_s     = search(data,r_find_s)
    mm_e     = search(data,r_find_e)
    str_s    = data[ (mm_s.stop+1) : (mm_e.start-1) ]
    #print(str_s)
    return str_s
end


function extractNamedSections(data::String, name::String)

    named_data::Dict{String,String} = Dict{String,String}()

    last_pos::Int64 = 1
    while(true)
        r_find_s = Regex("$(name)Start:{.+?}")
        ms = search(data,r_find_s,last_pos)
        if(ms.stop<0)
            break
        end
        last_pos = ms.stop+1
        # find section name..
        sss = data[ms]
        sss = replace(sss,"$(name)Start:","")
        sss = replace(sss,"{","")
        sss = replace(sss,"}","")
        sec_name = strip(sss)

        # find end..
        r_find_e = Regex("$(name)End:{$(sec_name)}")
        ms2 = search(data,r_find_e)

        str_s = data[ (ms.stop+1) : (ms2.start-1) ]
        #println(str_s)
        named_data[sec_name] = str_s
    end
    #println(named_data)
    return named_data
end

#removes comments, strips everything etc..
function preprocessLine(li::String)
    li = replace(li,"\r","")
    lis = split(li,"//")
    if(length(li)>1)
        li = String(lis[1])
    end
    li = strip(li)
    return li
end



function parseFlux(id::String, rs::String)
    sp1 = split(rs,"-->")
    rs  = sp1[1]
    ps  = sp1[2]

    rr  = split(rs,"+")
    pp  = split(ps,"+")

    rasd = Dict{String,RAS}()
    pasd = Dict{String,PAS}()

    # keep track of number of reactands of specific species:
    species_count_r::LC{String} = LC(String)
    for ri in rr
        (stoich::Float64,sid::String,mat::String) = parseReactand(String(ri))
        species_count_r = add( species_count_r , sid )
        if( (!isempty(mat)) && (stoich!=1) ); error("Non-one stoich. for reactand with atom transition map"); end
        # now parse atom map
        if(!isempty(mat))
            mats = split(mat,"-")
            sci::Int64  = 1
            for zs::String in mats
                sas1::SAS = SAS(sid,sci)
                #ras1::RAS = RAS( sas1 , Int64(round(species_count_r.map[sid])) )
                ras1::RAS = RAS( id , sas1 , Int64(round(species_count_r[sid])) )
                rasd[zs] = ras1
                sci = sci + 1
            end
        end
        #println("$(rasd)")
    end
    species_count_p::LC{String} = LC(String)
    for ri in pp
        (stoich::Float64,sid::String,mat::String) = parseReactand(String(ri))
        species_count_p = add( species_count_p , sid )
        if( (!isempty(mat)) && (stoich!=1) ); error("Non-one stoich. for reactand with atom transition map"); end
        # now parse atom map
        if(!isempty(mat))
            mats = split(mat,"-")
            sci::Int64  = 1
            for zs::String in mats
                sas1::SAS = SAS(sid,sci)
                #pas1::PAS = PAS( sas1 , Int64(round(species_count_p.map[sid])) )
                pas1::PAS = PAS( id , sas1 , Int64(round(species_count_p[sid])) )
                pasd[zs] = pas1
                sci = sci + 1
            end
        end
        #println("$(pasd)")
    end
    # Create atom transition map:
    tm::Dict{RAS,PAS} = Dict{RAS,PAS}()
    if( haskey(species_count_r,"*") || haskey(species_count_p,"*") )
        # no tm..
    else
        # compute tm
        if( Set(keys( rasd )) != Set(keys( pasd )) ); error("Parsing TM failed.. : $(id) : $(rs)"); end
        #keys_r = sort( keys( rasd ) )
        #keys_p = sort( keys( pasd ) )
        for ki in keys(rasd)
            tm[ rasd[ki] ] = pasd[ki]
        end
    end
    fx::Flux = Flux(id,species_count_r,species_count_p,tm)
    return fx
end

#returns string representation of flux
function write_atom_transition(f::Flux)
    transition_inv = map(reverse,f.transition)

    atom_placeholders = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" # i.e. exporting map with > 52 atoms will fail..
    str_ed::IOBuffer = IOBuffer()

    educts = sort(collect(keys(f.E)))
    map_atom_identifiers = Dict{RAS,String}()
    atom_counter = 1
    for ei in educts
        #find sas of ei:
        sas_ei = filter( s -> isequal( s.sas.species , ei ) , collect(keys(f.transition)) )
        #find max slot pos:
        num_slots = maximum( map( s -> s.sas.pos , sas_ei ) )

        for zs in 1:f.E[ei]
            print( str_ed , ei)
            print( str_ed , "(")
            for zx in 1:num_slots
                print( str_ed , atom_placeholders[atom_counter:atom_counter] )
                map_atom_identifiers[ RAS(f.id,SAS(ei,zx),zs) ] = atom_placeholders[atom_counter:atom_counter]
                atom_counter += 1
                if(zx<num_slots)
                    print( str_ed , "-")
                end
            end
            print( str_ed , ") + ")
        end
    end
    str_educts = String(take!(str_ed))
    str_educts = str_educts[1:end-2]

    products = sort(collect(keys(f.P)))
    str_prod::IOBuffer = IOBuffer()
    for pi in products
        #find sas of pi:
        sas_pi = filter( s -> isequal( s.sas.species , pi ) , collect(values(f.transition)) )
        #find max slot pos:
        num_slots = maximum( map( s -> s.sas.pos , sas_pi ) )

        for zs in 1:f.P[pi]
            print( str_prod , pi)
            print( str_prod , "(")
            for zx in 1:num_slots
                # find correct identifier:
                pas_i = PAS(f.id,SAS(pi,zx),zs)
                str_i = map_atom_identifiers[ transition_inv[pas_i] ]
                print( str_prod , str_i )
                if(zx<num_slots)
                    print( str_prod , "-")
                end
            end
            print( str_prod , ") + ")
        end
    end
    str_products = String(take!(str_prod))
    str_products = str_products[1:end-2]

    return string( str_educts , " --> " , str_products )
end

#return stoichiometry, "species id" and the atom transition map part
#i.e.:
# "sA()" --> (1.0,"sA","")
# " 4.0 sB " --> (4.0,"sB","")
# " sC(a-b-X) " --> (1.0,"sC","a-b-X")
function parseReactand(rs::String)
    debug = 0

    stoich::Float64 = 1.0
    ri1 = strip(rs)
    ri2 = split( rs , "(" )
    if(size(ri2,1)==1)
        ri1 = string(ri1,"()")
    end
    #extract species sid
    if(debug>0); println("ri: $(ri1)") ; end
    sid::String = ri1[ 1:(searchindex(ri1,"(")-1)]
    sid = strip(sid)

    # try to split sid..
    sid_s = split(sid,r"\s+")
    if(length(sid_s)>2); error("parseReactand: error for : \"$(rs)\""); end
    if(length(sid_s)==1); end
    if(length(sid_s)==2)
        sid = sid_s[2]
        stoich = float(String(sid_s[1]))
    end
    if(debug>0); println("stoich: $(stoich) sid: $(sid)") ; end
    #species_count = add( species_count , sid )
    #println(species_count)
    # extract atom map
    mat_a = match( r"\((.*)\)" , ri1 )
    mat = mat_a[1]
    mat = replace(mat,"(","")
    mat = replace(mat,")","")
    if(debug>0); println("mat: $(mat)"); end
    return (stoich,sid,mat)
end



function parseMeasurementSet(ms::String)

    mss = split(ms,"\n")

    in_mdset::Bool = false # are we currently parsing measurements

    parsed_measurements = Vector{MassDistributionMeasurement}()
    current_ms::MassDistributionMeasurement = MassDistributionMeasurement( Set{SAS}() , Dict{Int64,Array{Float64,1}}() )# = MassDistributionMeasurement(Set{SAS}(),Dict{}())

    for li in mss
        li = String(li)
        li = preprocessLine(li)
        if(isempty(li)); continue; end

        #println("Process line: $(String(li))")

        if(in_mdset)
            if(li=="}")
                push!(parsed_measurements,current_ms)
                in_mdset = false
                continue
            else
                lis = split(li,":")
                lis_a = lis[1]
                if(search(lis_a,"M=")[2]<1); error("Parsing of Mass wrong syntax.."); end
                lis_a = replace( lis_a , "M=" , "" )
                imass = parse(Int64,lis_a)
                # parse 3 floats for lis_b
                lis_b = lis[2]
                lis_b_s = split(lis_b,",")
                lis_b_values::Array{Float64,1} = [ parse(Float64,lis_b_s[1]) ; parse(Float64,lis_b_s[2]) ; parse(Float64,lis_b_s[3]) ]
                current_ms.data[imass] = lis_b_values
            end
        else
            lis = split(li,"{")
            lis = lis[1]
            liss = split(lis,"[")
            liss_a = liss[1]
            ms_name::String = strip(liss_a)
            liss_b = liss[2]
            liss_b = replace(liss_b,"]","")
            svals = split( liss_b ,"," )
            sass::Set{SAS} = Set{SAS}()
            for svi in svals
                push!( sass , SAS( ms_name , parse(Int64,svi) ) )
            end
            current_ms = MassDistributionMeasurement( sass , Dict{Int64,Array{Float64,1}}() )
            in_mdset = true
        end
    end
    return parsed_measurements
end


function parseSubstrateMixture(sm_str::String)

    sm_s = split(sm_str,"\n")

    parsed_species = Vector{ConfiguredSEMUMixture}()

    for li in sm_s
        li = String(li)
        li = preprocessLine(li)
        if(isempty(li)); continue; end

        #println("Process line: $(String(li))")

        li_s = split(li,":")
        if(length(li_s)!=2); error("Syntax Error in Substrate Mixture Def. : $(String(li))"); end
        li_a = strip(li_s[1]) # species id

        li_b = strip(li_s[2])
        li_b_s = split(li_b,";")

        csasm::ConfiguredSEMUMixture = ConfiguredSEMUMixture()

        for zi=1:length(li_b_s)
            li_b_s[zi] = strip(li_b_s[zi])
            #println(li_b_s[zi])
            mixci = split(li_b_s[zi],"=")
            if(length(mixci)!=2); error("Syntax Error in Substrate Mixture Def. : $(String(li))"); end
            mixci_a = mixci[1]
            mixci_a = replace(String(mixci_a),"(","")
            mixci_a = replace(String(mixci_a),")","")
            mixci_a_s = split(mixci_a,",")
            config_i::Array{Int64,1} = Array{Int64,1}()
            for isi in mixci_a_s
                push!(config_i,parse(Int64,isi))
            end
            sas_i::Array{SAS,1} = Array{SAS,1}()
            for zj=1:length(mixci_a_s)
                push!(sas_i,SAS(li_a,zj))
            end
            cr::ConfiguredSEMU = ConfiguredSEMU(Set{SAS}(sas_i),config_i)
            csasm[cr] = parse(Float64,mixci[2])
        end
        push!( parsed_species , csasm)
    end
    return parsed_species
end


function parseLinConstr(data::String)

    lcs = split(data,"\n")

    lin_constr = Vector{LinConstr}()

    for li in lcs
        li = String(li)
        li = preprocessLine(li)
        if(isempty(li)); continue; end
        #println("Process line: $(String(li))")
        push!( lin_constr , parse_linconstr(li) )
    end

    return lin_constr
end


function get_configured_reactands( fd::FluxiData , substr_conf::String)
    sc = fd.substrate_configs[substr_conf]
    csemus::Vector{ConfiguredSEMU} = Vector{ConfiguredSEMU}()
    for sci in sc
        append!( csemus , collect(keys(sci)) )
    end
    #print(csemus)

    all_sas::Vector{SAS} = Vector{SAS}()
    for cse in csemus
        append!( all_sas , cse.semu)
    end
    return Set(all_sas)
end
