

"""
FWSim_Julia contains the created julia code for the FW-Sim
"""
mutable struct FWSim_Julia <: FWSim
    code::String
    fluxes::Vector{String} # defines the order of fluxes which have to be passed to the function, is equal to the flux order in the supplied network
    semu_mass_ratios::Vector{ASTEMUMassRatio} # defines the order of defined semus
    target_mass_ratios::Vector{ASTEMUMassRatio} # defines the order of result semus , is equal to the supplied vector
end



""""
function ast_to_juliacode_01( data::FluxiData , ast::ASTEMUSim , target_smrs::Vector{ASTEMUMassRatio} ; fast_lu=true)

 compiles the ast to juliacode
 creates a julia_fw_sim object with:
-------------------------------------------------------------------------------
 .code has structure:
 "function( semu_mass_ratios::Vector{Float64} ,  v::Vector{Float64} )
     return md
   end"
-------------------------------------------------------------------------------
 the other fields define the order of the two supplied vectors and the result
 vector (see julia_fw_sim struct documentation for details)



 fast_lu : if true, then we use lu-decomp without pivoting. this SHOULD BE
           NO PROBLEM, IF supplied fluxes are NOT VERY CLOSE TO ZERO..
"""
function ast_to_juliacode_01( net::Network , ast::ASTEMUSim , target_smrs::Vector{ASTEMUMassRatio} ; fast_lu=true)
    flag_printDebugOutput = false

    code_header::IOBuffer = IOBuffer()

    # 1. create function "header"
    print(code_header , "function( semu_mass_ratios::Vector{Float64} ,  v::Vector{Float64} )\n\n")


    # here we add all that we find in the X matrix, in any step
    semu_mass_ratios_computed::Set{ASTEMUMassRatio} = Set{ASTEMUMassRatio}()

    # 3. step in, start with processing the sim steps..
    code_computation::IOBuffer = IOBuffer()

    req_constants::Dict{ASTNumNode,String} = Dict{ASTNumNode,String}()
    for siz in sort( collect( keys(ast.sim_steps) ) )
        ss_cnt = 0
        for ss in ast.sim_steps[siz]
            ss_cnt+=1
            print(code_computation, "#----Init A---------------------------------------------------------------------------------\n")
            print( code_computation , to_code_ast_matrix( "sA_$(siz)_$(ss_cnt)" , ss.A , req_constants ) )
            print(code_computation, "#----Init B---------------------------------------------------------------------------\n")
            print( code_computation , to_code_ast_matrix( "sB_$(siz)_$(ss_cnt)" , ss.B , req_constants ) )
            print(code_computation, "#----Init Y---------------------------------------------------------------------------\n")
            print( code_computation , to_code_ast_matrix( "sY_$(siz)_$(ss_cnt)" , ss.Y , req_constants ) )

            print(code_computation, "\n#------COMPUTE!-----------------------------------------------------------------------\n")
            # now add the computation of the matrix Y, and then read out the newly computed values
            if(!fast_lu)
                print( code_computation , "sX_$(siz)_$(ss_cnt) = sA_$(siz)_$(ss_cnt) \\ ( sB_$(siz)_$(ss_cnt) * sY_$(siz)_$(ss_cnt) )\n" )
            end
            if(fast_lu)
                print( code_computation , "( sA_L_$(siz)_$(ss_cnt) , sA_U_$(siz)_$(ss_cnt) ) = lu(sA_$(siz)_$(ss_cnt),Val{false})\n" )
                print( code_computation , "sX_pre_$(siz)_$(ss_cnt) = sA_L_$(siz)_$(ss_cnt) \\ ( sB_$(siz)_$(ss_cnt) * sY_$(siz)_$(ss_cnt) )\n" )
                print( code_computation , "sX_$(siz)_$(ss_cnt) = sA_U_$(siz)_$(ss_cnt) \\ sX_pre_$(siz)_$(ss_cnt) \n\n" )
            end

            if(flag_printDebugOutput)
                s_ma = "sA_$(siz)_$(ss_cnt)"
                s_mb = "sB_$(siz)_$(ss_cnt)"
                s_my = "sY_$(siz)_$(ss_cnt)"
                s_mx = "sX_$(siz)_$(ss_cnt)"

                println( code_computation , "print( \"STEP: A = \\n \") \n" )
                println( code_computation , "print($(s_ma))")
                println( code_computation , "print(\"\\n\")")
                println( code_computation , "print( \"STEP: B = \\n \") \n" )
                println( code_computation , "print($(s_mb))")
                println( code_computation , "print(\"\\n\")")
                println( code_computation , "print( \"STEP: Y = \\n \") \n")
                println( code_computation , "print($(s_my))")
                println( code_computation , "print(\"\\n\")")
                println( code_computation , "print( \"STEP: X = \\n \") \n")
                println( code_computation , "print($(s_mx))")
                println( code_computation , "print(\"\\n\")")
            end

            print(code_computation, "#-------------------------------------------------------------------------------------\n")
            # now read out the newly computed values from X:
            for zib = 1:size(ss.X.m,2)
                for zia = 1:size(ss.X.m,1)
                    smri = ss.X.m[zia,zib]
                    @assert typeof(smri) == ASTEMUMassRatio
                    var_name_si = create_var_name(smri)
                    req_constants[smri] = var_name_si
                    push!( semu_mass_ratios_computed , smri )
                    print(code_computation, "$(var_name_si) = sX_$(siz)_$(ss_cnt)[$(zia),$(zib)]\n")
                end
            end
            print(code_computation, "#-------------------------------------------------------------------------------------\n")
        end
    end

    # 4. now we add the required variable initializations
    # 4.1 collect the required mass distribution ratios:
    code_varinit::IOBuffer = IOBuffer()

    rc_flux = filter( x1 -> typeof(x1)==ASTFluxRate , collect(keys(req_constants)) )
    rc_emr  = filter( x1 -> typeof(x1)==ASTEMUMassRatio , collect(keys(req_constants)) )
    rc_flux = convert(Vector{ASTFluxRate},collect(rc_flux))
    rc_emr  = convert(Vector{ASTEMUMassRatio},collect(rc_emr))



    rc_emr_needed::Vector{ASTEMUMassRatio} = collect( setdiff( Set(rc_emr) , Set(semu_mass_ratios_computed) ) )
    rc_emr_needed                          = sort(rc_emr_needed)
    print("\n\nRequired EMRs: \n");
    print(code_varinit, "#-----Init EMRs-------------------------------------------------------------------------------------\n")
    emri_cnt = 0
    for emri in rc_emr_needed
        emri_cnt += 1
        println("$(emri)")
        print(code_varinit,"$(req_constants[emri]) = semu_mass_ratios[$(emri_cnt)]\n")
    end

    rc_flux = sort(rc_flux) # NOTE: this does not change anything, just makes the resulting code look nicer..
    print("\n\nRequired Fluxes: \n");
    print(code_varinit, "#-----Init Fluxes-----------------------------------------------------------------------------------\n")
    for fluxii in rc_flux
        println("$(fluxii)")
        #fluxpos = find( xfi -> xfi==fluxii.flux , data.net.fluxes )
        fluxpos = find( xfi -> xfi==fluxii.flux , net.fluxes )
        @assert length(fluxpos)==1
        print(code_varinit,"$(req_constants[fluxii]) = v[$(fluxpos[1])]\n")
    end
    print(code_varinit, "#---------------------------------------------------------------------------------------------------\n")

    # 5. collect the results
    code_results::IOBuffer = IOBuffer()
    print(code_results, "#---RESULTS-----------------------------------------------------------------------------------------\n")
    cnt_tsr = 0
    print(code_results , "target_smrs = Array{Float64,1}($(length(target_smrs)))\n")
    for t_semu in target_smrs
        cnt_tsr +=1
        print(code_results , "target_smrs[$(cnt_tsr)] = $(req_constants[t_semu])\n")
    end
    print(code_results, "#---------------------------------------------------------------------------------------------------\n")
    print(code_results,"return target_smrs\n")
    print(code_results, "#---------------------------------------------------------------------------------------------------\n")

    #str_total = string( takebuf_string(code_header) ,  takebuf_string(code_header) ,  takebuf_string(code_computation) )
    code_total = string( String(take!(code_header)) ,  String(take!(code_varinit)) ,  String(take!(code_computation) ) , String(take!(code_results) ) , "end" )

    #return FWSim_Julia(code_total,data.net.fluxes,rc_emr_needed,[])
    return FWSim_Julia(code_total,net.fluxes,rc_emr_needed,target_smrs)
end

function to_code_ast_matrix(matrix_name::String , m::ASTMatrix ,  req_constants::Dict{ASTNumNode,String} )
    sb::IOBuffer = IOBuffer()
    #print(sb,"$(matrix_name)::Array{Float64,2} = Array{Float64,2}($(size(m.m,1)),$(size(m.m,2))) \n")
    print(sb,"$(matrix_name)::Array{Float64,2} = zeros(Float64,$(size(m.m,1)),$(size(m.m,2))) \n")

    # fill the matrix..
    for zj=1:size(m.m,2)
        for zi=1:size(m.m,1)
            # evaluate cell..
            if(isassigned(m.m,zi,zj))
                str_cell = to_code_ast_num(m.m[zi,zj],req_constants)
                print(sb,"$(matrix_name)[$(zi),$(zj)] = $(str_cell) \n")
            end
        end
    end
    return String( take!( sb ) )
end

function to_code_ast_num( n::ASTSymmOp , req_constants::Dict{ASTNumNode,String})
    str = "( "
    op_str = "$(n.op)"
    for zi=1:length(n.args)
        str = string( str , to_code_ast_num(n.args[zi] , req_constants) , " " , op_str )
    end
    str = str[1:length(str)-1]
    str = string(str," ) ")
    return str
end

function to_code_ast_num( n::ASTRealValue , req_constants::Dict{ASTNumNode,String} )
    return " $(n.val) "
end

function to_code_ast_num( n::ASTFluxRate , req_constants::Dict{ASTNumNode,String} )
    if( ! haskey( req_constants , n ) )
        req_constants[n] = create_var_name(n)
    end
    return " $(req_constants[n]) "
end
function to_code_ast_num( n::ASTEMUMassRatio , req_constants::Dict{ASTNumNode,String} )
    if( ! haskey( req_constants , n ) )
        req_constants[n] = create_var_name(n)
    end
    return " $(req_constants[n]) "
end

function create_var_name( n::ASTFluxRate)
    return string("v_",n.flux)
end

function create_var_name( n::ASTEMUMassRatio)
    @assert isinsinglespecies(n.semu)
    ss::String = string("mr_",collect(n.semu)[1].species,"_")
    n_semu_pos = sort( collect( map( x1 -> x1.pos , n.semu ) ) )
    for seip in n_semu_pos
        ss = string( ss , seip , "_" )
    end
    #ss = ss[1:length(ss)-1]
    mass_i =  round(Int , n.mass )
    ss = string( ss , "M_$(mass_i)")
    return String(ss)
end






function benchmark_matrix_solving_01()

    # create lots of toy problems to solve..
    # i.e. diagonal-heavy matrix to invert

    m_sizes        = [100;20;6]
    m_density      = 0.25
    m_diagheavy_a  = 0.1
    m_diagheavy_b  = 1.5
    # A: ms[1] x ms[1] , B: ,s[1] x ms[2] , Y: s[2] x x[3]
    # diagheavy: to the diagonal of m we add
    # dh_a and dh_b times the abs value of the row

    n_tests = 1000

    test_A = Vector{Array{Float64,2}}(n_tests)
    test_B = Vector{Array{Float64,2}}(n_tests)
    test_Y = Vector{Array{Float64,2}}(n_tests)

    test_X_A = Vector{Array{Float64,2}}(n_tests)
    test_X_B = Vector{Array{Float64,2}}(n_tests)
    test_X_C = Vector{Array{Float64,2}}(n_tests)
    test_X_D = Vector{Array{Float64,2}}(n_tests)

    failed_problems = Vector{Array{Array{Float64,2}}}()

    # prepare..
    println("\nPrepare Test:\n")
    for zi=1:n_tests
        Ai = rand(m_sizes[1],m_sizes[1])
        Ai[Ai.<(1.-m_density)] = 0
        for zx=1:size(Ai,1); Ai[zx,zx] = Ai[zx,zx] + m_diagheavy_a  + m_diagheavy_b*sum(abs.(Ai[zx,:])); end
        Bi = rand(m_sizes[1],m_sizes[2])
        Bi[Bi.<(1.-m_density)] = 0
        Yi = rand(m_sizes[2],m_sizes[3])
        Yi[Yi.<(1.-m_density)] = 0

        Ai = Ai * 20 * rand()
        Bi = Bi * 20 * rand()
        Yi = Yi * 20 * rand()

        test_A[zi] = Ai
        test_B[zi] = Bi
        test_Y[zi] = Yi
        if(rem(zi,n_tests/20)==0); print('.'); end
    end

    println("\nRun Test:\n")

    if(false)
        print("Test A: ")
        tic()
        for zi=1:n_tests
            test_X_A[zi] = inv(test_A[zi]) * test_B[zi] * test_Y[zi]
            if(rem(zi,n_tests/20)==0); print('.'); end
        end
        println("")
        toc()
        println("DONE!")

        print("Test B: ")
        tic()
        for zi=1:n_tests
            test_X_A[zi] = test_A[zi] \ ( test_B[zi] * test_Y[zi] )
            if(rem(zi,n_tests/20)==0); print('.'); end
        end
        println("")
        toc()
        println("DONE!")
    end
    if(true)
        print("Test B: ")
        tic()
        for zi=1:n_tests
            test_X_B[zi] = test_A[zi] \ ( test_B[zi] * test_Y[zi] )
            if(rem(zi,n_tests/20)==0); print('.'); end
        end
        println("\n")
        toc()
        println("DONE!")
    end
    if(true)
        print("Test C: ")
        tic()
        for zi=1:n_tests
            (A_L,A_U,A_P) = lu(test_A[zi])
            x = Array{Float64,2}(1,1)

            x_pre = A_L \ ( ( test_B[zi] * test_Y[zi] )[A_P,:] )
            x     = A_U \ x_pre
            test_X_C[zi] = x
            #test_X_B[zi] = test_A[zi] \ ( test_B[zi] * test_Y[zi] )

            if(rem(zi,n_tests/20)==0); print('.'); end
        end
        println("")
        toc()
        println("DONE!")
    end
    if(true)
        print("\nTest D: ")
        tic()
        for zi=1:n_tests
            (A_L,A_U,A_P) = lu(test_A[zi],Val{false})
            x_pre = A_L \ ( test_B[zi] * test_Y[zi] )
            x     = A_U \ x_pre
            test_X_D[zi] = x
            #test_X_B[zi] = test_A[zi] \ ( test_B[zi] * test_Y[zi] )

            if(rem(zi,n_tests/20)==0); print('.'); end
        end
        println("")
        toc()
        println("DONE!")
    end

    # compute max error in D:
    f_dist = (x,y) -> sum( ((x.-y).^2)[:] )
    errors_d = Array{Float64,1}(length(test_X_B))
    for zi=1:length(test_X_B)
        errors_d[zi] = f_dist(test_X_B[zi],test_X_D[zi])
    end
    println("Max Error in D: $(maximum(errors_d))")
end


function test_lu_solving_01()
    a1 = rand(1000,1000)
    b1 = rand(1000,12)

    (a_l,a_u,a_p) = lu(a1)
    x_pre = a_l \ b1[a_p,:]; x = a_u \ x_pre;

    (a2_l,a2_u) = lu(a1,Val{false})
    x_pre_2 = a2_l \ b1; x_2 = a2_u \ x_pre_2;

    print("\n\n")
    #println(a1\b1)
    #println(x)
    #println(x_2)
    x_ref = a1\b1
    f_dist = (x,y) -> sum( ((x.-y).^2)[:] )

    println("Diff 1: $(f_dist(x_ref,x))")
    println("Diff 2: $(f_dist(x_ref,x_2))")
end
