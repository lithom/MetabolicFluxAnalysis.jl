

struct FTBLConstr
    kind::Char # can be F or C , (D are skipped)
    value::Float64
end

struct FTBL_MDMeasurement_Entry
    semu::SpeciesEMU
    mass::Int64
    data::Vector{Float64}
end

function FTBL_MDMeasurement_Entry(semu::SpeciesEMU,mass::Int64,mean::Float64,stdev::Float64)
    return FTBL_MDMeasurement_Entry(semu,mass,[mean;stdev])
end

struct FTBL_Data
    fluxes::Vector{Flux}
    constraints_net::Dict{String,FTBLConstr}
    constraints_xch::Dict{String,FTBLConstr}

    eq_constr_net::Vector{LinConstr}
    eq_constr_xch::Vector{LinConstr}
    ineq_constr_net::Vector{LinConstr}
    ineq_constr_xch::Vector{LinConstr}

    substrate_conf::Dict{ConfiguredSEMU,Float64}
    md_measurements::Vector{FTBL_MDMeasurement_Entry}

    flux_measurements::Dict{String,Vector{Float64}}
end

struct FTBL_Import
    fluxidata::FluxiData
    wx_transform::WX_Transform
    wx_constraints_fluxes::Vector{LinConstr}
    wx_constraints_eq_and_ineq::Vector{LinConstr}
    wx_constraints_flux_measurements::Vector{LinConstr}
    wx_constraints_fluxes_proposed::Vector{LinConstr} # contains the "F" recommended values vas lin. constraints
end

"""
  ftbl_to_fluxi(data::FTBL_Data)

This function creates a fluxi dataset including the fluxi network from the parsed FTBL data.
It also creates a WX-Parametrization for the network, on which the constraints from the FTBL_Data are applied.

# NOTES

1. For fluxes with non-zero xch flux, we create a fw and a bw flux (name of forward flux is original name, backward flux has _bw postfix)

2. Create In/Out-Flows: To "open" the network, we create Influxes for:
    (i)  the metabolites in "substrate_conf"
  and we create Outfluxes for:
    (i)  all "dead-end" metabolites, i.e. metabolites with only in-fluxes     (I guess this is what 13C flux does?!)

"""
function ftbl_to_fluxi(data::FTBL_Data)

    fluxes = copy(data.fluxes)

    # create the fluxi network..


    # 1.  Collect all species:
    # 1.1 collect sas of all metabolites:
    species = Set{Species}()
    all_sas = reduce( (x,y) -> union(x,y) , map( x -> union( map( ri->ri.sas , keys(x.transition)) , map(pi->pi.sas ,values(x.transition)))  ,  data.fluxes ) )
    # 1.2 create species:
    all_species_ids = unique( map( sas -> sas.species , all_sas) )
    for si in all_species_ids
        sas_si = filter( sas -> sas.species==si , all_sas )
        num_slots = maximum( map( sas -> sas.pos , sas_si ) )
        push!(species,Species(si,num_slots))
    end

    # 2. Check for dead-end species and create out-fluxes for those..
    all_ras = reduce( (x,y) -> union(x,y) , map( x -> keys(x.transition) , data.fluxes ) ) # !! Do this on the original data.fluxes
    # if a species never appears in a ras, then it is a dead-end metabolite..
    all_ras_species = unique( map( ri -> ri.sas.species , all_ras ) )
    dead_end_metabolites = setdiff( map( si -> si.id , species) , all_ras_species )
    print("dead end metabolites: \n")
    for dmi in dead_end_metabolites;
        print("$(dmi) -> create Ouflux..\n")
        # create outflux:
        push!( fluxes , MetabolicFluxAnalysis.parseFlux(string("Outflux_",dmi) , string( dmi , " --> * " ) ) )
    end




    # 2.b collect all influx sas:
    influx_sas = reduce( (x,y)->union(x,y) , sas.(collect(keys(data.substrate_conf))) )
    #print(influx_sas)


    # 3. Loop over all fluxes and check which need bw flux:
    # NOTE: also fluxes which lead into dead-end metabolites do NOT have a bw flux (I guess?!)
    # NOTE: also the in-fluxes do NOT have a bw flux (I guess?!)
    # NOTE: when deciding, that a flux has no bw flux, based on the "C" constraint, then
    # we HAVE TO REMOVE the constraint, else it will fail..
    for fi in data.fluxes
        needs_bw_flux = true
        if( haskey(data.constraints_xch,fi.id) )
            if( data.constraints_xch[fi.id].kind=='C' )
                if(data.constraints_xch[fi.id].value == 0)
                    needs_bw_flux = false
                end
            end
        end

        # check for fluxes to dead end metabolites:
        if( length( setdiff( dead_end_metabolites , map( si->si.sas.species ,values(fi.transition))) ) < length(dead_end_metabolites) )
            print("found flux to dead-end metabolite -> no reverse flux : $(fi.id)\n")
            needs_bw_flux = false
        end
        if(  any( in.( map( ri -> ri.sas , collect( keys(fi.transition ) ) ) , [  influx_sas  ] ) ) )
            print("found input flux -> no reverse flux : $(fi.id)\n")
            needs_bw_flux = false
        end
        if(needs_bw_flux)
            print("create reverse flux for $(fi.id)\n")
            insert!( fluxes  , find( f -> f==fi , fluxes )[1]+1 ,  create_reverse_flux(fi) )
        end
    end



    # 4a. Prepare substrates:
    # just collect the different semus
    # !! AND ADD THE INFLUX !!
    substr_semus = unique( map( sci -> sci.semu , collect(keys(data.substrate_conf)) ) )
    sorted_substr = Vector{ConfiguredSEMUMixture}()
    for ssemu in substr_semus
        csemu = Dict{ConfiguredSEMU,Float64}()
        for sci in filter( sci -> sci.semu==ssemu , collect(keys(data.substrate_conf)) )
            csemu[sci] = data.substrate_conf[sci]
        end
        push!( sorted_substr , csemu )
    end
    substr_configs = Dict{String,Vector{ConfiguredSEMUMixture}}()
    substr_configs["ftbl_default"] = sorted_substr
    # 4b. add influxes:
    influx_species = unique( reduce( (x,y) -> union(x,y) , map( si -> map(s->s.species , si.semu) , collect(keys(data.substrate_conf))  ) ) )
    for si in influx_species
        print("$(si) -> create Influx..\n")
        push!( fluxes , parseFlux( string("Influx_",si) , string(" * --> ",si) ))
    end



    # 5. The network should be complete now, so create it:
    species_db = Dict{String,Species}()
    for si in species; species_db[si.id] = si; end
    fluxes_db  = Dict{String,Flux}()
    for fi in fluxes; fluxes_db[fi.id] = fi; end

    net = Network( map(f -> f.id,fluxes) , map(f -> f.id,collect(species)) , fluxes_db , species_db )



    # 6. Create mass distribution measurements:
    # again, just collect the different semus
    mdm_semus   = unique( map( mdi -> mdi.semu , data.md_measurements ) )
    sorted_mdms = Vector{MassDistributionMeasurement}()
    for mdsi in mdm_semus
        ddi = Dict{Int64,Vector{Float64}}()
        for mmi in filter( mi -> mi.semu==mdsi , data.md_measurements)
            ddi[mmi.mass] = mmi.data
        end
        push!( sorted_mdms , MassDistributionMeasurement(mdsi,ddi) )
    end
    mdm_configs = Dict{String,Vector{MassDistributionMeasurement}}()
    mdm_configs["ftbl_default"] = sorted_mdms


    # 7. WX Parametrization and constraints..
    (p_pwx,wxt)  = MetabolicFluxAnalysis.parametrization_wiechert_exchange( net , Array{MetabolicFluxAnalysis.LinConstr,1}() , 200. )

    # 8. init the WX constraints
    # We can just create normal LinConstr, with "WX-Fluxes", i.e. with [x] [n] notation:
    constraints_c          = Vector{LinConstr}()
    constraints_proposed   = Vector{LinConstr}() # contains the "F" values..

    # 8.1 process the data.constraints_net
    for fi in keys(data.constraints_net)
        if(data.constraints_net[fi].kind == 'C')
            push!( constraints_c , parse_linconstr( string( fi,"[n] = ", data.constraints_net[fi].value) ) )
        else
            push!( constraints_proposed , parse_linconstr( string( fi,"[n] = ", data.constraints_net[fi].value) ) )
        end
    end
    # 8.2 process the data.constraints_xch
    for fi in keys(data.constraints_xch)
        if(is_wx_flux(net,wxt,fi)) # NO CONSTRAINT,IF CONSTR. TO ZERO, THEN WE DID NOT ADD THE XCH FLUX!!
            if( data.constraints_xch[fi].kind == 'C' )
                push!( constraints_c , parse_linconstr( string( fi,"[x] = ", data.constraints_xch[fi].value) ) )
            else
                push!( constraints_proposed , parse_linconstr( string( fi,"[x] = ", data.constraints_xch[fi].value) ) )
            end
        else
            # don't add anything
        end
    end
    # 8.3 process the data.eq_constraints_net and data.ineq_constriants_net
    # we just change the fids in the lin constraint and add the "[n]"
    net_constraints = [data.eq_constr_net;data.ineq_constr_net]
    wx_constraints  = Vector{LinConstr}()
    for nci in net_constraints
        new_rhs = Dict{String,Float64}()
        for rhi in keys(nci.rhs)
            if(rhi=="1"); new_rhs["1"]=nci.rhs[rhi];else; new_rhs[string(rhi,"[n]")]=nci.rhs[rhi];end
        end
        new_lhs = Dict{String,Float64}()
        for lhi in keys(nci.lhs)
            if(lhi=="1"); new_lhs["1"]=nci.lhs[lhi];else; new_lhs[string(lhi,"[n]")]=nci.lhs[lhi];end
        end
        push!( wx_constraints , LinConstr(nci.constr,new_lhs,new_rhs) )
    end
    # 8.4 process the data.eq_constraints_xch and data.ineq_constriants_xch
    # we just change the fids in the lin constraint and add the "[x]"
    # NOTE: We have to pay attention for [x] fluxes that we did not include, be-
    # cause they are constrained to zero. Therefore, we check for all LC entries,
    # if they are wx-parametrized, and if not, we replace the [x] flux with zero.
    xch_constraints = [data.eq_constr_xch;data.ineq_constr_xch]
    for xci in xch_constraints
        new_rhs = Dict{String,Float64}()
        for rhi in keys(xci.rhs)
            if(rhi=="1"); new_rhs["1"]=xci.rhs[rhi]
            else
                if(is_wx_flux(net,wxt,rhi))
                    new_rhs[string(rhi,"[x]")]=xci.rhs[rhi]
                else
                    # nothing to do, is just zero
                end
            end
        end
        new_lhs = Dict{String,Float64}()
        for lhi in keys(xci.lhs)
            if(lhi=="1")
                new_lhs["1"]=xci.lhs[lhi]
            else
                if(is_wx_flux(lhi))
                    new_lhs[string(lhi,"[x]")]=xci.lhs[lhi]
                else
                    # nothing to do, is just zero
                end
            end
        end
        push!( wx_constraints , LinConstr(xci.constr,new_lhs,new_rhs) )
    end

    # 8.5 add the constraints from .flux_measurements
    # for simplicity, we just add them as equality constraints (we could also derive ineq. constr easily, but..)
    # NOTE: as these flux measurements are usually for in/out fluxes (which do not have WX-parametrization), we
    # add the flux_id without [x/n] postfix..
    # We may have to fix this if we have flux measurements for WX-Fluxes..
    flux_measurements_wx_constraints = Vector{LinConstr}()
    for ki in keys( data.flux_measurements )
        push!( flux_measurements_wx_constraints , parse_linconstr( string( ki , " = " , data.flux_measurements[ki][1] ) ) )
    end

    # Ok, now we have everything to create the fluxi dataset:
    fdata = FluxiData(net,mdm_configs,substr_configs,Dict{String,Vector{LinConstr}}() )


    wx_constraints = [wx_constraints;]
    return FTBL_Import(fdata,wxt,constraints_c, wx_constraints,flux_measurements_wx_constraints,constraints_proposed)
    #return (net,substr_configs,mdm_configs,wxt,MetabolicFluxAnalysis.create_all_wx_fluxnames(net,wxt) , constraints_c , constraints_proposed , wx_constraints )
end

"""
  parse_ftbl_data(data::IOStream)

  parses the data from .ftbl data and returns FTBL_Data object

  #NOTES
  This function parses all sections, except: PROJECT , LABELING_MEASUREMENTS , PEAK_MEASUREMENTS , OPTIONS

"""
function parse_ftbl_data(data::IOStream)

    # loop over lines:
    lines = readlines(data)

    line_pos = 1

    #
    # parsing status can be:
    #
    # 10: in NETWORK section, waiting for new flux
    # 11: in NETWORK section, waiting for atom transition map
    #
    # 20: in FLUXES section, before NET or XCH
    # 21: in FLUXES NET section
    # 22: in FLUXES XCH section
    # 30: in EQUALITIES section, before NET or XCH
    # 31: in EQUALITIES NET section
    # 32: in EQUALITIES XCH section
    # 41: in INEQUALITIES section, before NET or XCH
    # 41: in INEQUALITIES NET section
    # 42: in INEQUALITIES XCH section
    # 50: in MASS_SPECTROMETRY section
    # 60: in LABEL_INPUT section
    # 70: in FLUX_MEASUREMENTS section
    parsing_status = -1
    parsing_done   = false


    # data while constructing fluxes:
    fluxes             = Vector{Flux}()
    mfluxes            = Vector{String} # holds the parsed tsv of the first part of a flux

    # data while collecting constraints
    constraints_net = Dict{String,FTBLConstr}()
    constraints_xch = Dict{String,FTBLConstr}()

    eq_constr_net   = Vector{LinConstr}()
    eq_constr_xch   = Vector{LinConstr}()
    ineq_constr_net   = Vector{LinConstr}()
    ineq_constr_xch   = Vector{LinConstr}()

    flux_measurements = Dict{String,Vector{Float64}}() # 2-element arrays, with [value ; deviation]

    # data for collecting label_input data:
    # we just put all lines that we find in here, and post-process it to create the substrate mixture
    configured_semus = Dict{ConfiguredSEMU,Float64}()
    current_config_metabolite = ""

    # data for collecting mass distribution measurements:
    md_measurements = Vector{FTBL_MDMeasurement_Entry}()
    current_md_semu = SpeciesEMU()


    while( line_pos <= length(lines) )
        li_a = lines[line_pos]

        li = preprocess_line(li_a)
        if( test_empty_line(li) )
            line_pos += 1
            continue
        end

        if(startswith(li,"NETWORK"))
            parsing_status = 10
            line_pos       = line_pos+2 # skip that bloody "info" line
            continue
        end
        if(startswith(li,"FLUXES"))
            parsing_status = 20
            line_pos       = line_pos+1
            continue
        end
        if(startswith(li,"EQUALITIES"))
            parsing_status = 30
            line_pos       = line_pos+1
            continue
        end
        if(startswith(li,"INEQUALITIES"))
            parsing_status = 40
            line_pos       = line_pos+1
            continue
        end
        if(startswith(li,"MASS_SPECTROMETRY"))
            parsing_status = 50
            line_pos       = line_pos+2 # skip that bloody "info" line
            continue
        end
        if(startswith(li,"LABEL_INPUT"))
            parsing_status = 60
            line_pos       = line_pos+2 # skip that bloody "info" line
            continue
        end
        if(startswith(li,"FLUX_MEASUREMENTS"))
            parsing_status = 70
            line_pos       = line_pos+2 # skip that bloody "info" line
            continue
        end
        # and parse the sections that we do not parse:
        if(startswith(li,"LABEL_MEASUREMENTS"))
            parsing_status = -1
            line_pos       = line_pos+2 # skip that bloody "info" line
            continue
        end
        if(startswith(li,"PEAK_MEASUREMENTS"))
            parsing_status = -1
            line_pos       = line_pos+2  # skip that bloody "info" line
            continue
        end
        if(startswith(li,"OPTIONS"))
            parsing_status = -1
            line_pos       = line_pos+2  # skip that bloody "info" line
            continue
        end




        if( parsing_status==20 || parsing_status==21 || parsing_status==22 )
            if( strip(li)=="NET"  )
                parsing_status = 21
                line_pos += 2 # skip that bloody "info" line
                continue
            end
            if( strip(li)=="XCH"  )
                parsing_status = 22
                line_pos += 2 # skip that bloody "info" line
                continue
            end
        end
        if( parsing_status==30 || parsing_status==31 || parsing_status==32 )
            if( strip(li)=="NET"  )
                parsing_status = 31
                line_pos += 2 # skip that bloody "info" line
                continue
            end
            if( strip(li)=="XCH"  )
                parsing_status = 32
                line_pos += 2 # skip that bloody "info" line
                continue
            end
        end
        if( parsing_status==40 || parsing_status==41 || parsing_status==42 )
            if( strip(li)=="NET"  )
                parsing_status = 41
                line_pos += 2 # skip that bloody "info" line
                continue
            end
            if( strip(li)=="XCH"  )
                parsing_status = 42
                line_pos += 2 # skip that bloody "info" line
                continue
            end
        end




        # now add the parse code..
        if(parsing_status==10)
            #mflux = matchall(r"\w+",li)
            mfluxes = parse_tsv_line(li)
            parsing_status = 11
            line_pos += 1
            continue;
        end
        if(parsing_status==11)
            tmap = parse_tsv_line(li)
            # now use fluxi parser to construct flux.. ;)
            fluxi_flux = ""
            if(!isempty(mfluxes[3]))
                fluxi_flux = string( fluxi_flux , mfluxes[3] , "(" , conv_tmap(tmap[3]) , ")" )
            end
            if(!isempty(mfluxes[4]))
                if(!isempty(mfluxes[3])); fluxi_flux = string(fluxi_flux , " + "); end
                fluxi_flux = string( fluxi_flux , mfluxes[4] , "(" , conv_tmap(tmap[4]) , ")" )
            end
            fluxi_flux = string(fluxi_flux , " --> " )
            if(!isempty(mfluxes[5]))
                fluxi_flux = string( fluxi_flux , mfluxes[5] , "(" , conv_tmap(tmap[5]) , ")" )
            end
            if(length(mfluxes)>=6)
                if(!isempty(mfluxes[6]))
                    if(!isempty(mfluxes[5])); fluxi_flux = string(fluxi_flux , " + "); end
                    fluxi_flux = string( fluxi_flux , mfluxes[6] , "(" , conv_tmap(tmap[6]) , ")" )
                end
            end
            push!( fluxes , MetabolicFluxAnalysis.parseFlux( mfluxes[2] , fluxi_flux) )

            parsing_status = 10
            line_pos += 1
            continue
        end

        if(parsing_status==21)
            fct = parse_tsv_line(li)
            if( fct[4]=="F" || fct[4]=="C" )
                constraints_net[ fct[3] ] = FTBLConstr( fct[4][1] , parse(Float64,fct[5]) )
            elseif( fct[4]=="D" )
            else
                print("Skipped constraint with unclear type: $(li)\n")
            end
            line_pos += 1
            continue
        end
        if(parsing_status==22)
            #NOTE: xch constraints to zero we DO include, because based on those
            #we later decide, how the wx transform looks..
            fct = parse_tsv_line(li)
            if(length(fct)<4)
                print("sumtingwong: $(li)\n")
            end
            if( fct[4]=="F" || fct[4]=="C" )
                #if( parse(Float64,fct[5])==0 && fct[4]=="C")
                    #print("Zero xch constraint sorted out\n")
                #else
                    constraints_xch[ fct[3] ] = FTBLConstr( fct[4][1] , parse(Float64,fct[5]) )
                #end
            elseif( fct[4]=="D" )
            else
                print("Skipped constraint with unclear type: $(li)")
            end
            line_pos += 1
            continue
        end

        if(parsing_status==31)
            cpt = parse_tsv_line(li)
            push!( eq_constr_net , MetabolicFluxAnalysis.parse_linconstr( string(cpt[3], " = ",cpt[4]) ) )
            line_pos += 1
            continue
        end
        if(parsing_status==32)
            cpt = parse_tsv_line(li)
            push!( eq_constr_xch , MetabolicFluxAnalysis.parse_linconstr( string(cpt[3], " = ",cpt[4]) ) )
            line_pos += 1
            continue
        end
        if(parsing_status==41)
            #NOTE!!! The ftbl. notation is weird. Constraint is in OPPOSITE direction of what you read!!!
            cpt = parse_tsv_line(li)
            #if( cpt[4]=="<=" || cpt[4]=="<" )
            if( cpt[4]==">=" || cpt[4]==">" )
                push!( ineq_constr_xch , MetabolicFluxAnalysis.parse_linconstr( string(cpt[3], " < ",cpt[5]) ) )
            #elseif( cpt[4]==">=" || cpt[4]==">" )
        elseif( cpt[4]=="<=" || cpt[4]=="<" )
                push!( ineq_constr_xch , MetabolicFluxAnalysis.parse_linconstr( string(cpt[3], " > ",cpt[5]) ) )
            else
                print("ineq constraint with unclear comp. $(li)")
            end
            line_pos += 1
            continue
        end
        if(parsing_status==42)
            cpt = parse_tsv_line(li)
            if( cpt[4]=="<=" || cpt[4]=="<" )
                push!( ineq_constr_xch , MetabolicFluxAnalysis.parse_linconstr( string(cpt[3], " < ",cpt[5]) ) )
            elseif( cpt[4]==">=" || cpt[4]==">" )
                push!( ineq_constr_xch , MetabolicFluxAnalysis.parse_linconstr( string(cpt[3], " > ",cpt[5]) ) )
            else
                print("ineq constraint with unclear comp. $(li)")
            end
            line_pos += 1
            continue
        end

        if(parsing_status==70)
            fmt = parse_tsv_line(li)
            flux_measurements[fmt[2]] = [ parse(Float64,fmt[3]) ; parse(Float64,fmt[4]) ]
            line_pos += 1
            continue
        end

        if(parsing_status==60)
            llt = parse_tsv_line(li)
            if(!isempty(llt[2]))
                current_config_metabolite = llt[2]
            end
            configured_semus[ ftbl_create_configured_semu(current_config_metabolite,llt[3]) ] = parse(Float64,llt[4])
            line_pos += 1
            continue
        end

        if(parsing_status==50)
            mdt = parse_tsv_line(li)
            if(!isempty(mdt[2]))
                current_md_semu = constructSpeciesEMU( mdt[2] , map( x -> parse(Int64,x) , split(mdt[3],",") ) )
            end
            push!( md_measurements , FTBL_MDMeasurement_Entry( current_md_semu , parse(Int64,mdt[4]) , parse(Float64,mdt[5]) , parse(Float64,mdt[6]) ) )
            line_pos += 1
            continue
        end

        # if we arrive here, nothing parsed this line..
        print("unparsed line: $(li)\n")
        line_pos += 1
        continue
    end

    #return ( fluxes,constraints_net,constraints_xch , (eq_constr_net,eq_constr_xch,ineq_constr_net,ineq_constr_xch) , configured_semus , md_measurements )
    return FTBL_Data(fluxes,constraints_net,constraints_xch , eq_constr_net,eq_constr_xch,ineq_constr_net,ineq_constr_xch , configured_semus , md_measurements , flux_measurements )
end


function preprocess_line(a::String)
    m_comment = match(r"//", a)
    if(isa(m_comment,Void))
        return a
    else
        a         = a[1:m_comment.offset-1]
        return a
    end
end

function test_empty_line(a::String)
    a2 = a[:]
    m_has_word_characters = match(r"\w",a2)
    return isa(m_has_word_characters,Void)
end


function parse_tsv_line(li::String)

    ti = find( x -> x=='\t',li)
    all_indeces = [1;ti;length(li)+1]

    tabbed::Vector{String} = Vector{String}(length(all_indeces)-1)

    for zi in 1:(length(all_indeces)-1)
        tabbed[zi] = li[ all_indeces[zi]+1 : all_indeces[zi+1]-1]
    end
    return tabbed
end
# convergts #ABC into A-B-C
function conv_tmap(orig::String)
    tmap = ""
    for zi=2:length(orig)
        tmap = string(tmap,orig[zi])
        if(zi<length(orig)); tmap = string(tmap,"-"); end
    end
    return tmap
end

# Converts "sA" , "#001100"
# into the corresponding ConfiguredSEMU
function ftbl_create_configured_semu( met::String , conf_str::String)
    conf_str = strip(conf_str)
    csemu    = constructSpeciesEMU( met , collect( 1:length(conf_str)-1 ) )
    conf     = map( x -> parse(Int,x) , collect(conf_str[2:end]))
    return ConfiguredSEMU(csemu,conf)
end
