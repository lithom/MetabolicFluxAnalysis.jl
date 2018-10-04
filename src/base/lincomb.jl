
import Base.show


#struct LC{T}
#    map::Dict{T,Float64}
#end


LC{T} = Dict{T,Float64}

function LC(dt::DataType)
     #return LC(Dict{dt,Float64}())
     return LC{dt}()
end

function parseLC(lcs::String)
    debug=0
    #lc = LC{String}(Dict{String,Float64}())
    lc = LC{String}()
    ss = split(lcs,"+")
    for si in ss
        sis= strip(si)
        #sx = split( sis , r"\s+" )
        sx = split( sis , r"[\s\*]+" )
        vx::String = ""
        nx::Float64 = 1.
        if(size(sx,1)==2)
            if(debug>0); println("$(sis) --> $(sx[1]) $(sx[2])") ; end
            nx = float(sx[1])
            vx = strip(sx[2])
        else
            if(debug>0); println("$(sis) --> $(sx[1])") ; end
            vx = strip(sx[1])
        end
        #lci = LC{String}( Dict(vx => nx) )
        lci = LC{String}( vx => nx )
        #lci.map[vx] = nx
        lc = add(lc,lci)
    end
    return lc
end

function add{T}(a::LC{T},b::LC{T})
    #keys_a::Array{T,1} = collect(keys(a.map))
    #keys_b::Array{T,1} = collect(keys(b.map))
    keys_a::Array{T,1} = collect(keys(a))
    keys_b::Array{T,1} = collect(keys(b))

    keys_ab::Array{T,1} = unique([keys_a;keys_b])
    #map_ab  = Dict{T,Float64}()
    map_ab  = LC{T}()
    for x in keys_ab
        #vi::Float64 = ( haskey(a.map,x) ) ? a.map[x] : 0
        #vi         += ( haskey(b.map,x) ) ? b.map[x] : 0
        vi::Float64 = ( haskey(a,x) ) ? a[x] : 0
        vi         += ( haskey(b,x) ) ? b[x] : 0
        map_ab[x] = vi
    end
    #return LC{T}(map_ab)
    return map_ab
end

function add{T}(a::LC{T},b::T)
    return add(a,b,1.0)
end

function add{T}(a::LC{T},b::T,count::Float64)
    lcb = LC{T}( b => count )
    return add(a,lcb)
end

function mult{T}(a::LC{T},b::Float64)
    lc_new = LC{T}()
    for ki in keys(a)
        lc_new[ki] = b * a[ki]
    end
    return lc_new
end

"""
LinConstr has fields constr::Char, lhs::LC{String} and rhs::LC{String}
and represents a linear constraint
.constr: can be "=", "<" , ">"
convention: "1" in lhs/rhs is the key for the non-variable-multiplied value
"""
struct LinConstr
    constr::Char
    lhs::LC{String}
    rhs::LC{String}
end

function Base.show(io::IO,lc::LC{String})
    cnt=0
    for ki in keys(lc)
        cnt+=1
        #if(lc[ki]==0)
            # skip entry
        #else
            #print( io , @sprintf("%5.3d",lc[ki]) )
            print( io , lc[ki])
            print( io , string( " ",ki ) )
            if(cnt<length(lc)); print(io," + "); end
        #end
    end
end
function Base.show(io::IO,lc::LinConstr)
    show(io,lc.lhs)
    print(io,string(" ",lc.constr," "))
    show(io,lc.rhs)
end

function resolve_lc_s(n::Int64,d_s::Dict{String,Int64},lc::LC{String})
    lcrow::Array{Float64,1} = zeros(n)
    for ei in keys(lc)
        if(!haskey(d_s,ei)); error("LC::resolve_lc : error : no element $(ei)"); end
        lcrow[d_s[ei]] = lc[ei]
    end
    return lcrow
end

function resolve_lc{T}(n::Int64,d_s::Dict{T,Int64},lc::LC{T})
    lcrow::Array{Float64,1} = zeros(n)
    for ei in keys(lc)
        if(!haskey(d_s,ei)); error("LC::resolve_lc : error : no element $(ei)"); end
        lcrow[d_s[ei]] = lc[ei]
    end
    return lcrow
end

"""
parses linear constraints
comparison operators are: '=' , '<' , '>'
on left and right hand side there must be strings which can be parse
by parseLC
"""
function parse_linconstr(str)
    iseq = contains(str,"=")
    isle = contains(str,"<")
    isge = contains(str,">")

    n_elg = (iseq?1:0)+(isle?1:0)+(isge?1:0)
    if(n_elg!=1); error(string("Unparsable constraint: ",str)); end

    op::Char = '?'

    if(iseq) ; op ='=' ; end
    if(isle) ; op ='<' ; end
    if(isge) ; op ='>' ; end

    str_ab = split(str,string(op))
    a = String( str_ab[1] )
    b = String( str_ab[2] )

    #lhs = parseLC(a)
    #rhs = parseLC(b)
    lhs = parse_lc_expr(strip(a))
    rhs = parse_lc_expr(strip(b))

    return LinConstr(op,lhs,rhs)
end

"""
parse_lc_expr(str::String)
returns LC{String}

parses formulas which are sums like:
" 5 va + 4 bxy - 2 + 1 va "
" 3 "
" 4 v "
" 3 * x + 2 * abc + 3"
"""
function parse_lc_expr(str::String)
    # enry "1" maps to the "not-multiplied-with-var" part
    val = LC{String}("1" => 0)

    str = strip(str);
    if(isempty(str)); return val; end

    if(contains(str,"+"))
        split_pos = search(str,"+")[1]
        va = parse_lc_expr( String( str[1:split_pos-1] ) )
        vb = parse_lc_expr( String( str[split_pos+1:length(str)] ) )
        return add(va,vb)
    else
        if(contains(str,"-"))
            split_pos = search(str,"-")[1]
            va = parse_lc_expr( String( str[1:split_pos-1] ) )
            vb = parse_lc_expr( String( str[split_pos+1:length(str)] ) )
            return add(va, mult( vb , -1.0 ) )
        end
    end

    ## now down to parsing literals..
    #str_s = split(str,"\s+")
    str_s = split(str,r"[\s\*]+")
    if(length(str_s)==1)
        str_s = strip(string(str_s[1]))
        #str_s = String(str_s)
        #rmatch = ismatch(r"\d*\.?\d+",str_s)
        parsed_float = false
        try
            pf = parse(Float64,str_s)
            val["1"] = pf
            parsed_float = true
        end
        if(!parsed_float)
            val[str_s] = 1.0
        end
        return val
    else
        if(length(str_s)>2); error("Error parsing constraint: "+str); end
        v_num = parse(Float64,str_s[1])
        v_var = strip(String(str_s[2]))
        val[v_var] = v_num
        return val
    end
end
