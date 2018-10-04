
using OpenCL


struct OCLContext
    ocl_device::OpenCL.cl.Device
    ocl_ctx::OpenCL.cl.Context
    ocl_queue::OpenCL.cl.CmdQueue
end

function OCLContext()
    device, ctx, queue = OpenCL.cl.create_compute_context()
    return OCLContext(device, ctx, queue)
end

function OCLContext(device::OpenCL.cl.Device)
    ctx   = cl.Context(device);
    queue = cl.CmdQueue(ctx)
    return OCLContext(device,ctx,queue)
end

"""
FWSim_OpenCL contains the created OpenCL code for the FW-Sim

Structure of the code is a single kernel with 3 parameters:
__kernel void fwsim(__global const float *semus,
                  __global const float *fluxes,
                  __global float *mds_out)

semus[ i*length(semu_mass_ratios) ] is the first element of the semus data
for sample i.
fluxes[ i*length(fluxes) ] is the first element of the flux rate data
for sample i.
"""
mutable struct FWSim_OpenCL <: FWSim
    code::String
    fluxes::Vector{String} # defines the order of fluxes which have to be passed to the function, is equal to the flux order in the supplied network
    semu_mass_ratios::Vector{ASTEMUMassRatio} # defines the order of defined semus
    target_mass_ratios::Vector{ASTEMUMassRatio} # defines the order of result semus , is equal to the supplied vector

    required_matrix_size_x::Int
    required_matrix_size_y::Int

    ocl_context::Nullable{OCLContext}
    k::Nullable{OpenCL.cl.Kernel}
    p::Nullable{OpenCL.cl.Program}
end
function FWSim_OpenCL( code::String, fluxes::Vector{String}, semu_mass_ratios::Vector{ASTEMUMassRatio} , target_mass_ratios::Vector{ASTEMUMassRatio} , required_matrix_size_x::Int , required_matrix_size_y::Int  )
    return FWSim_OpenCL(code,fluxes,semu_mass_ratios,target_mass_ratios, required_matrix_size_x,required_matrix_size_y,Nullable{OCLContext}() , Nullable{OpenCL.cl.Kernel}() , Nullable{OpenCL.cl.Program}() )
end

function ast_to_opencl_01( net::Network , ast::ASTEMUSim , target_smrs::Vector{ASTEMUMassRatio} )

    flag_include_debug::Bool = false

    # 1. create function "header"
    code_header::IOBuffer = IOBuffer()
    print(code_header ,"     __kernel void fwsim(__global const float *semus,
                          __global const float *fluxes,
                          __global float *mds_out)
        {
          const int gid = get_group_id(0) + get_group_id(1);
          const int idx = get_local_id(0);
          const int idy = get_local_id(1);

          //printf(\"\\n group_id=%d idx=%d idy=%d glob_x=%d glob_y=%d \\n\",gid,idx,idy,get_global_id(0),get_global_id(1));

          // Here we store the sizes of the given matrices in a step
          int mAm = 0;
          int mAn = 0;
          int mBm = 0;
          int mBn = 0;
          int mXm = 0;
          int mXn = 0;
          int mYm = 0;
          int mYn = 0;


        \n\n")

    # 1.b Init the local matrices needed for the single steps..
    req_size_A = maximum( map( x -> maximum( collect(map( xi -> prod(size(xi.A.m)) , x))) , values(ast.sim_steps) ) )
    req_size_B = maximum( map( x -> maximum( collect(map( xi -> prod(size(xi.B.m)) , x))) , values(ast.sim_steps) ) )
    req_size_X = maximum( map( x -> maximum( collect(map( xi -> prod(size(xi.X.m)) , x))) , values(ast.sim_steps) ) )
    req_size_Y = maximum( map( x -> maximum( collect(map( xi -> prod(size(xi.Y.m)) , x))) , values(ast.sim_steps) ) )
    req_size_BY = maximum( map( x -> maximum( collect(map( xi -> size(xi.B.m,1) * size(xi.Y.m,2) , x))) , values(ast.sim_steps) ) )
    # 1.c Create local matrices:
    print(code_header , "__local float mA[$(req_size_A)]; \n")
    print(code_header , "__local float mA_L[$(req_size_A)]; \n") # lu decomp, l part
    print(code_header , "__local float mB[$(req_size_B)]; \n")
    print(code_header , "__local float mX[$(req_size_X)]; \n")
    print(code_header , "__local float mY[$(req_size_Y)]; \n")
    print(code_header , "__local float mBY[$(req_size_BY)]; \n")

    # here we add all that we find in the X matrix, in any step
    semu_mass_ratios_computed::Set{ASTEMUMassRatio} = Set{ASTEMUMassRatio}()

    # 3. step in, start with processing the sim steps..
    code_computation::IOBuffer = IOBuffer()

    req_constants::Dict{ASTNumNode,String} = Dict{ASTNumNode,String}()
    for siz in sort( collect( keys(ast.sim_steps) ) )
        ss_cnt = 0
        for ss in ast.sim_steps[siz]
            ss_cnt+=1
            print(code_computation, "//----Init A---------------------------------------------------------------------------------\n")
            print(code_computation, "if(idx==0 && idy==0){\n")
            print( code_computation , to_opencl_code_ast_matrix( "mA" , ss.A , req_constants ) )
            print(code_computation, "}\n")
            print(code_computation, "//----Init B---------------------------------------------------------------------------\n")
            print(code_computation, "if(idx==1 && idy==0){\n")
            print( code_computation , to_opencl_code_ast_matrix( "mB" , ss.B , req_constants ) )
            print(code_computation, "}\n")
            print(code_computation, "//----Init Y---------------------------------------------------------------------------\n")
            print(code_computation, "if(idx==2 && idy==0){\n")
            print( code_computation , to_opencl_code_ast_matrix( "mY" , ss.Y , req_constants ) )
            print(code_computation, "}\n")

            # !! We have to set the matrix size indicators on work item level!!
            print(code_computation,"mAm = $(size(ss.A.m,1)); \n")
            print(code_computation,"mAn = $(size(ss.A.m,2)); \n")
            print(code_computation,"mBm = $(size(ss.B.m,1)); \n")
            print(code_computation,"mBn = $(size(ss.B.m,2)); \n")
            print(code_computation,"mYm = $(size(ss.Y.m,1)); \n")
            print(code_computation,"mYn = $(size(ss.Y.m,2)); \n")



            # might not be needed..
            print(code_computation, "barrier(CLK_LOCAL_MEM_FENCE);\n")

            if(flag_include_debug) # some debug code.. :)
                print( code_computation ,  "\nif(idx==0&&idy==0){printf(\"\\nBefore Computation:\");}\n")
                print( code_computation ,  ocl_print_matrix( "mA" , size(ss.A.m,1) , size(ss.A.m,2)) )
                print( code_computation ,  ocl_print_matrix( "mB" , size(ss.B.m,1) , size(ss.B.m,2)) )
                print( code_computation ,  ocl_print_matrix( "mY" , size(ss.Y.m,1) , size(ss.Y.m,2)) )
            end

            print(code_computation, "\n//------COMPUTE!-----------------------------------------------------------------------\n")

            #print(code_computation, to_opencl_code_ludecomp_2("mA",size(ss.A.m,1)) )
            print(code_computation, to_opencl_code_ludecomp_2_1d("mA",size(ss.A.m,1)) )
            print(code_computation, "//----------------------------------\n")
            #print(code_computation, to_opencl_code_matrixmul("mBY","mB","mY", ss.B , ss.Y ) )
            print(code_computation, to_opencl_code_matrixmul_1d("mBY","mB","mY", ss.B , ss.Y ) )
            print(code_computation,"barrier(CLK_LOCAL_MEM_FENCE);\n")

            if(flag_include_debug) # some debug code.. :)
                print( code_computation ,  "\nif(idx==0&&idy==0){printf(\"\\nResult LU Decomp:\");}\n")
                print( code_computation ,  ocl_print_matrix( "mA" , size(ss.A.m,1) , size(ss.A.m,2)) )
                print( code_computation ,  ocl_print_matrix( "mA_L" , size(ss.A.m,1) , size(ss.A.m,2)) )
            end

            if(flag_include_debug) # some debug code.. :)
                print( code_computation ,  "\nif(idx==0&&idy==0){printf(\"\\nafter Multiply B and Y:\");}\n")
                print( code_computation ,  ocl_print_matrix( "mBY" , size(ss.B.m,1) , size(ss.Y.m,2)) )
            end

            print(code_computation, to_opencl_code_solve_L_system( "mA_L" , "mBY" , size(ss.A.m,1) , size(ss.A.m,2) , size(ss.Y.m,2) ))
            print(code_computation,"barrier(CLK_LOCAL_MEM_FENCE);\n")

            if(flag_include_debug) # some debug code.. :)
                print( code_computation ,  "\nif(idx==0&&idy==0){printf(\"\\nafter L Solve:\");}\n")
                print( code_computation ,  ocl_print_matrix( "mBY" , size(ss.B.m,1) , size(ss.Y.m,2)) )
            end

            print(code_computation, to_opencl_code_solve_U_system( "mA" , "mBY" , size(ss.A.m,1) , size(ss.A.m,2) , size(ss.Y.m,2) ))
            print(code_computation,"barrier(CLK_LOCAL_MEM_FENCE);\n")

            if(flag_include_debug) # some debug code.. :)
                print( code_computation ,  "\nif(idx==0&&idy==0){printf(\"\\nafter U Solve:\");}\n")
                print( code_computation ,  ocl_print_matrix( "mBY" , size(ss.B.m,1) , size(ss.Y.m,2)) )
            end


            #    print( code_computation , "( sA_L_$(siz)_$(ss_cnt) , sA_U_$(siz)_$(ss_cnt) ) = lu(sA_$(siz)_$(ss_cnt),Val{false})\n" )
            #    print( code_computation , "sX_pre_$(siz)_$(ss_cnt) = sA_L_$(siz)_$(ss_cnt) \\ ( sB_$(siz)_$(ss_cnt) * sY_$(siz)_$(ss_cnt) )\n" )
            #    print( code_computation , "sX_$(siz)_$(ss_cnt) = sA_U_$(siz)_$(ss_cnt) \\ sX_pre_$(siz)_$(ss_cnt) \n\n" )

            print(code_computation, "//-------------------------------------------------------------------------------------\n")
            #print(code_computation, "if(idx==0&&idy==0){ \n") # !! ALL WORK ITEMS HAVE TO DO THIS!!
            # now read out the newly computed values from X:
            for zib = 1:size(ss.X.m,2)
                for zia = 1:size(ss.X.m,1)
                    smri = ss.X.m[zia,zib]
                    @assert typeof(smri) == ASTEMUMassRatio
                    var_name_si = create_var_name(smri)
                    req_constants[smri] = var_name_si
                    push!( semu_mass_ratios_computed , smri )
                    #print(code_computation, "$(var_name_si) = sX_$(siz)_$(ss_cnt)[$(zia),$(zib)]\n")
                    print(code_computation, "$(var_name_si) = mBY[ $( (zia-1) + (zib-1) * size(ss.X.m,1)) ];\n") # !! -1
                end
            end
            #print(code_computation, "} \n")
            print(code_computation, "barrier(CLK_LOCAL_MEM_FENCE); \n")
            print(code_computation, "//-------------------------------------------------------------------------------------\n")
        end
    end

    # 4. now we add the required variable initializations
    # 4.1 collect the required mass distribution ratios:
    # code_varinit also sets the offsets for flux/semu/mratio data
    code_varinit::IOBuffer = IOBuffer()

    rc_flux = filter( x1 -> typeof(x1)==ASTFluxRate , collect(keys(req_constants)) )
    rc_emr  = filter( x1 -> typeof(x1)==ASTEMUMassRatio , collect(keys(req_constants)) )
    rc_flux = convert(Vector{ASTFluxRate},collect(rc_flux))
    rc_emr  = convert(Vector{ASTEMUMassRatio},collect(rc_emr))



    rc_emr_needed::Vector{ASTEMUMassRatio} = collect( setdiff( Set(rc_emr) , Set(semu_mass_ratios_computed) ) )
    rc_emr_needed                          = sort(rc_emr_needed)
    print("\n\nRequired EMRs: $(length(rc_emr_needed))  Required Fluxes: $(length(rc_flux))  Target MDs: $(length(target_smrs)) \n");
    print(code_varinit,"\n\n//Required EMRs: $(length(rc_emr_needed))  Required Fluxes: $(length(rc_flux))  Target MDs: $(length(target_smrs)) \n");
    print(code_varinit,"int offset_emrs   = 0; \n"); # CHANGE THIS IF EVERY WORK ITEM HAS DIFF. Substrate Config..
    print(code_varinit,"int offset_fluxes = gid * $(length(rc_flux)); \n");
    print(code_varinit,"int offset_target = gid * $(length(target_smrs)); \n");

    #print(code_varinit,"\n\nRequired EMRs: \n");
    print(code_varinit, "//-----Init EMRs-------------------------------------------------------------------------------------\n")
    emri_cnt = 0
    for emri in rc_emr_needed
        emri_cnt += 1
        println("$(emri)")
        print(code_varinit,"float $(req_constants[emri]) = semus[offset_emrs + $(emri_cnt-1)];\n") #!!!! -1
    end

    rc_flux = sort(rc_flux) # NOTE: this does not change anything, just makes the resulting code look nicer..
    #print("\n\nRequired Fluxes: \n");
    print(code_varinit, "//-----Init Fluxes-----------------------------------------------------------------------------------\n")
    for fluxii in rc_flux
        println("$(fluxii)")
        #fluxpos = find( xfi -> xfi==fluxii.flux , data.net.fluxes )
        fluxpos = find( xfi -> xfi==fluxii.flux , net.fluxes )
        @assert length(fluxpos)==1
        print(code_varinit,"float $(req_constants[fluxii]) = fluxes[ offset_fluxes + $(fluxpos[1]-1)];\n") # !!! -1
    end
    print(code_varinit, "//---------------------------------------------------------------------------------------------------\n")
    print(code_varinit, "//------Declare All Temporary EMR Float values--------------------------------------------------------------\n")
    all_emrs_without_initial_ones = setdiff( collect(keys(req_constants)) , union( rc_emr_needed , rc_flux ) )
    for emi in all_emrs_without_initial_ones
        print(code_varinit,"float $(req_constants[emi]);\n")
    end
    print(code_varinit, "//---------------------------------------------------------------------------------------------------\n")
    # 5. collect the results
    code_results::IOBuffer = IOBuffer()
    print(code_results, "//---RESULTS-----------------------------------------------------------------------------------------\n")
    print(code_results, "if(idx==0&&idy==0){ \n")
    cnt_tsr = 0
    #print(code_results , "target_smrs = Array{Float64,1}($(length(target_smrs)))\n")
    for t_semu in target_smrs
        cnt_tsr +=1
        print(code_results , "mds_out[ offset_target + $(cnt_tsr-1) ] = $(req_constants[t_semu]);\n") #!!!! -1
    end
    print(code_results, "}\n")
    print(code_results, "//---------------------------------------------------------------------------------------------------\n")
    #print(code_results,"return target_smrs\n")
    #print(code_results, "#---------------------------------------------------------------------------------------------------\n")

    #str_total = string( takebuf_string(code_header) ,  takebuf_string(code_header) ,  takebuf_string(code_computation) )
    code_total = string( String(take!(code_header)) ,  String(take!(code_varinit)) ,  String(take!(code_computation) ) , String(take!(code_results) ) , "}" )

    #return FWSim_Julia(code_total,data.net.fluxes,rc_emr_needed,[])
    #return code_total

    # find required matrix size..
    req_matrix_size = maximum( map( x -> maximum( collect(map( xi -> maximum( [ size(xi.A.m,1) ; size(xi.A.m,2) ;size(xi.B.m,1) ;size(xi.B.m,2) ;size(xi.Y.m,1) ;size(xi.Y.m,2) ] ) , x))) , values(ast.sim_steps) ) )

    return FWSim_OpenCL(code_total,net.fluxes,rc_emr_needed,target_smrs , req_matrix_size , req_matrix_size )
end


function ocl_mpos(m::ASTMatrix,x::Int,y::Int)
    return x*size(m.m,1)+y
end
function ocl_mpos(s::Int,x::Int,y::Int)
    return x*size(m.m,1)+y
end

function to_opencl_code_ast_matrix(matrix_name::String , m::ASTMatrix ,  req_constants::Dict{ASTNumNode,String} )
    sb::IOBuffer = IOBuffer()
    #print(sb,"$(matrix_name)::Array{Float64,2} = Array{Float64,2}($(size(m.m,1)),$(size(m.m,2))) \n")

    # this has to happen on work-item level!!
    #print(sb,"$(matrix_name)m = $(size(m.m,1)); \n")
    #print(sb,"$(matrix_name)n = $(size(m.m,2)); \n")

    # fill the matrix..
    for zj=0:size(m.m,2)-1
        for zi=0:size(m.m,1)-1
            # evaluate cell..
            if(isassigned(m.m,zi+1,zj+1))
                str_cell = to_opencl_code_ast_num(m.m[zi+1,zj+1],req_constants)
                print(sb,"$(matrix_name)[$(ocl_mpos(m,zj,zi))] = $(str_cell); \n")
            else
                print(sb,"$(matrix_name)[$(ocl_mpos(m,zj,zi))] = 0.0f; \n")
            end
        end
    end
    return String( take!( sb ) )
end

function to_opencl_code_ast_num( n::ASTSymmOp , req_constants::Dict{ASTNumNode,String})
    str = "( "
    op_str = "$(n.op)"
    for zi=1:length(n.args)
        str = string( str , to_opencl_code_ast_num(n.args[zi] , req_constants) , " " , op_str )
    end
    str = str[1:length(str)-1]
    str = string(str," ) ")
    return str
end

function to_opencl_code_ast_num( n::ASTRealValue , req_constants::Dict{ASTNumNode,String} )
    return " $(n.val)f "
end

function to_opencl_code_ast_num( n::ASTFluxRate , req_constants::Dict{ASTNumNode,String} )
    if( ! haskey( req_constants , n ) )
        req_constants[n] = create_var_name(n)
    end
    return " $(req_constants[n]) "
end
function to_opencl_code_ast_num( n::ASTEMUMassRatio , req_constants::Dict{ASTNumNode,String} )
    if( ! haskey( req_constants , n ) )
        req_constants[n] = create_var_name(n)
    end
    return " $(req_constants[n]) "
end


# the resulting code will compute lu decomposition of the given matrix.
# U will be in place of the original matrix, and L
#
# function assumes that we have a 2d group, at least as large as the matrix
#
function to_opencl_code_ludecomp( matrix_name::String , size::Int)
    sb::IOBuffer = IOBuffer()
    mn  = matrix_name;
    mnl = string(matrix_name,"_L")

    # 1. Init eye matrix for _L
    print(sb,"//LU-DECOMPOSITION (in place): $(matrix_name)_L * $(matrix_name) = $(matrix_name)\n")
    print(sb,"if(idx==0){
            for(int i=0; i<$(size*size); i++){ $(mnl)[i]=0f; }
            for(int i=0; i<$(size); i++){ $(mnl)[i*$(size)+i]=1f; }
        }\n")
    print(sb,"barrier(CLK_LOCAL_MEM_FENCE);\n")
    # 2. start with lu decomp:
    print(sb, "if(idx==0){\n")
    print(sb, "  for(int k=0;k<$(size);k++){\n")
    print(sb, "    for(int ii=k+1;ii<$(size);ii++){ $(mnl)[ k*$(size)+ii ] = $(mn)[ $(mnl)[ k*$(size)+ii ] / $(mn)[ k*$(size)+k ] ; } \n")

    print(sb, "    for(int ii=k+1;ii<$(size);ii++){ for(int ax=0;ax<$(size);ax++) { $(mn)[ ax*$(size)+ii ] -= $(mnl)[ k*$(size)+ii ] * $(mn)[ ax*$(size)+k ] ; } } \n")

    print(sb, "} }\n")

    return String( take!( sb ) )
end


# the resulting code will compute lu decomposition of the given matrix.
# U will be in place of the original matrix, and L
#
# function assumes that we have a 2d group, at least as large as the matrix
#
function to_opencl_code_ludecomp_2( matrix_name::String , size::Int)
    sb::IOBuffer = IOBuffer()
    mn  = matrix_name;
    mnl = string(matrix_name,"_L")
    s   = size

    # 1. Init eye matrix for _L
    print(sb,"//LU-DECOMPOSITION (in place): $(matrix_name)_L * $(matrix_name) = $(matrix_name)\n")
    print(sb,"if( idx<$(s) && idy<$(s) ){ $(mnl)[ idx+$(s)*idy ] = (idx==idy)?1.0f:0.0f; }\nbarrier(CLK_LOCAL_MEM_FENCE);\n")

    # 2. Loop over steps of lu decomp..
    for k in 0:s-1
        print(sb,"//LU-DECOMPOSITION STEP $(k) / $(s-1)\n")
        # process the L submatrix  (k+1:s) x (k)
        # process the U submatrix  (k+1:s) x (0:s)   //(k:s-1)


        # first handle L:
        # work item (0,idy) processes element (idy,k)
        print(sb, "if( idx==0 && idy > $(k) && idy < $(s) ){ \n");
        print(sb, "  $(mnl)[ idy + $(k*s) ] = $(mn)[ idy + $(k*s) ] / $(mn)[ $( k+k*s ) ]; \n");    # k*s is ok here
        print(sb, "} \n");
        print(sb, "barrier(CLK_LOCAL_MEM_FENCE);\n")

        # then handle U:
        # work item (idx,idy) processes element (idx,idy)
        print(sb, "if( idx>$(k) && idx < $(s) && idy < $(s) ){ \n");
        print(sb, "  $(mn)[ idx + $(s)*idy ] -= $(mnl)[ idx + $(s*k) ] * $(mn)[$(k)+ idy*$(s) ]; \n");
        print(sb, "} \n");
        print(sb, "barrier(CLK_LOCAL_MEM_FENCE);\n")

        if(false) # add some debug output..
            print(sb , "\nif(idx==0&&idy==0){ printf(\"After LU It:\\n\"); }\n")
            print(sb , ocl_print_matrix( mn , size , size))
            print(sb , ocl_print_matrix( mnl , size , size))
        end
    end
    print(sb,"//LU-DECOMPOSITION DONE!\n")
    return String( take!( sb ) )
end

# the resulting code will compute lu decomposition of the given matrix.
# U will be in place of the original matrix, and L
#
# function assumes that we have a 1d group, at least as large as the height of the matrix
#
function to_opencl_code_ludecomp_2_1d( matrix_name::String , size::Int)
    sb::IOBuffer = IOBuffer()
    mn  = matrix_name;
    mnl = string(matrix_name,"_L")
    s   = size

    # 1. Init eye matrix for _L
    print(sb,"//LU-DECOMPOSITION (in place): $(matrix_name)_L * $(matrix_name) = $(matrix_name)\n")
    #print(sb,"if( idx<$(s) && idy<$(s) ){ $(mnl)[ idx+$(s)*idy ] = (idx==idy)?1.0f:0.0f; }\nbarrier(CLK_LOCAL_MEM_FENCE);\n")
    print(sb,"if( idx<$(s)){\n")
    for idy_v in 0:s-1
         print(sb,"$(mnl)[ idx+$(s)*$(idy_v) ] = (idx==$(idy_v))?1.0f:0.0f; \n")
     end
     print(sb, "}\nbarrier(CLK_LOCAL_MEM_FENCE);\n" )


    # 2. Loop over steps of lu decomp..
    for k in 0:s-1
        print(sb,"//LU-DECOMPOSITION STEP $(k) / $(s-1)\n")
        # process the L submatrix  (k+1:s) x (k)
        # process the U submatrix  (k+1:s) x (0:s)   //(k:s-1)


        # first handle L:
        # work item (idx) processes element (idx,k)
        #print(sb, "if( idx==0 && idy > $(k) && idy < $(s) ){ \n");
        print(sb, "if( idx > $(k) && idx < $(s) && idy==0 ){ \n");
        print(sb, "  $(mnl)[ idx + $(k*s) ] = $(mn)[ idx + $(k*s) ] / $(mn)[ $( k+k*s ) ]; \n");    # k*s is ok here
        print(sb, "} \n");
        print(sb, "barrier(CLK_LOCAL_MEM_FENCE);\n")

        # then handle U:
        # work item (idx,idy) processes element (idx,idy)

        #print(sb, "if( idx>$(k) && idx < $(s) && idy < $(s) ){ \n");
        print(sb, "if( idx>$(k) && idx < $(s) && idy==0 ){ \n");
        for idy_v in 0:s-1
            #print(sb, "  $(mn)[ idx + $(s)*idy ] -= $(mnl)[ idx + $(s*k) ] * $(mn)[$(k)+ idy*$(s) ]; \n");
            print(sb, "  $(mn)[ idx + $(s)*$(idy_v) ] -= $(mnl)[ idx + $(s*k) ] * $(mn)[$(k)+ $(idy_v)*$(s) ]; \n");
        end
        print(sb, "} \n");
        print(sb, "barrier(CLK_LOCAL_MEM_FENCE);\n")

        if(false) # add some debug output..
            print(sb , "\nif(idx==0&&idy==0){ printf(\"After LU It:\\n\"); }\n")
            print(sb , ocl_print_matrix( mn , size , size))
            print(sb , ocl_print_matrix( mnl , size , size))
        end
    end
    print(sb,"//LU-DECOMPOSITION DONE!\n")
    return String( take!( sb ) )
end


# function assumes that we have a 2d group, at least as large as (mbm x myn)
#
function to_opencl_code_matrixmul( mby::String , mb::String , my::String , b::ASTMatrix , y::ASTMatrix )
    sb::IOBuffer = IOBuffer()
    mbm = size(b.m,1)
    myn = size(y.m,2)
    mbn = size(b.m,2)
    mym = size(y.m,1)

    # 1.
    print(sb,"//MATRIX MUL: $(mby) = $(mb) * $(my)   \n")
    print(sb,"//Work item (idx,idy) compute $(mby)[idx,idy]   \n")

    print(sb,"if( idx < $(mbm) && idy < $(myn) ) {  \n")
    # now this work item just computes this matrix element..
    # unroll.. ;)
    print(sb,"float mulsum = 0.0f;\n")
    for ni in 0:mbn-1
        print(sb,"mulsum += $(mb)[ $(mbm)*$(ni) + idx ] * $(my)[ idy * $(mym) + $(ni) ];\n")
    end
    print(sb,"$(mby)[idx+$(mbm)*idy] = mulsum;\n}\n")
    print(sb,"//MATRIX MUL COMPLETE! \n")

    return String( take!( sb ) )
end

# function assumes that we have a 1d group, at least as large as the height of b
#
function to_opencl_code_matrixmul_1d( mby::String , mb::String , my::String , b::ASTMatrix , y::ASTMatrix )
    sb::IOBuffer = IOBuffer()
    mbm = size(b.m,1)
    myn = size(y.m,2)
    mbn = size(b.m,2)
    mym = size(y.m,1)

    # 1.
    print(sb,"//MATRIX MUL: $(mby) = $(mb) * $(my)   \n")
    print(sb,"//Work item (idx) compute $(mby)[idx,:]   \n")

    #print(sb,"if( idx < $(mbm) && idy < $(myn) ) {  \n")
    print(sb,"if( idx < $(mbm) && idy < $(myn) ) {  \n")
    print(sb,"float mulsum = 0.0f;\n")
    #print(sb,"for(int idy_v=0; idy_v < $(myn);idy_v++ ) {  \n")
    for idy_v in 0:myn-1

        # now this iteration just computes this matrix element..
        # unroll.. ;)
        print(sb,"mulsum = 0.0f;\n")
        for ni in 0:mbn-1
            print(sb,"mulsum += $(mb)[ $(mbm)*$(ni) + idx ] * $(my)[ $(idy_v) * $(mym) + $(ni) ];\n")
        end
        print(sb,"$(mby)[idx+$(mbm)*$(idy_v)] = mulsum;\n")
    end
    print(sb,"}\n")

    print(sb,"//MATRIX MUL COMPLETE! \n")

    return String( take!( sb ) )
end

#
# Solves mL * x = mB
# solution will be in mB
#
# Not very parallel, work item (j) processes column j of mB
function to_opencl_code_solve_L_system( ml::String , mb::String ,  lm::Int  , ln::Int ,  bn::Int )
    sb::IOBuffer = IOBuffer()

    print(sb,"// Solve L System\n")
    #unroll..
    for i in 0:lm-1
        print(sb,"// Solve L System - Row $(i)\n")
        #print(sb,"if(idx==0 && idy<$(bn) ){\n")
        print(sb,"if(idx<$(bn) && idy==0 ){\n")
        print(sb,"  float sum_i = 0.0f;\n")
        for j in 0:(i-1)
            #print(sb,"  sum_i += $(mb)[ $(j) + $(lm)*idy ];\n")       # bm = lm
            #print(sb,"  sum_i += $(mb)[ $(j) + $(lm)*idy ] * $(ml)[ $(i) + $(j*lm) ] ;\n")       # bm = lm
            print(sb,"  sum_i += $(mb)[ $(j) + $(lm)*idx ] * $(ml)[ $(i) + $(j*lm) ] ;\n")       # bm = lm
        end
        #print(sb,"  $(mb)[ $(i) + $(lm)*idy ] = ($(mb)[ $(i) + $(lm)*idy ] - sum_i ) / $(ml)[ $(i+lm*i) ] ;\n")         # bm = lm
        print(sb,"  $(mb)[ $(i) + $(lm)*idx ] = ($(mb)[ $(i) + $(lm)*idx ] - sum_i ) / $(ml)[ $(i+lm*i) ] ;\n")         # bm = lm
        print(sb , "}\n")
        print(sb, "barrier(CLK_LOCAL_MEM_FENCE);\n")
    end
    print(sb,"// Solve L System DONE!\n")
    return String( take!( sb ) )
end

#
# Solves mU * x = mB
# In place, solution will be in mB
#
# Work item (j) processes column j of mB
function to_opencl_code_solve_U_system( mu::String , mb::String ,  um::Int  , un::Int ,  bn::Int )
    sb::IOBuffer = IOBuffer()

    print(sb,"// Solve U System\n")
    #unroll..
    for i_t in 0:(um-1) #(um-1):-1:0
        i = (um-1)-i_t
        print(sb,"// Solve U System - Row $(i)\n")
        #print(sb,"if(idx==0 && idy<$(bn) ){\n")
        print(sb,"if(idx<$(bn) && idy==0 ){\n")
        print(sb,"  float sum_i = 0.0f;\n")
        for j in ((un-1)-i_t+1):(un-1)  # un:-1:( (un-1)-i+1)  # 0:(i-1)
            #print(sb,"  sum_i += $(mb)[ $(j) + $(um)*idy ] * $(mu)[ $(i) + $(j*um) ] ;\n")       # bm = lm
            print(sb,"  sum_i += $(mb)[ $(j) + $(um)*idx ] * $(mu)[ $(i) + $(j*um) ] ;\n")       # bm = lm
        end
        #print(sb,"  $(mb)[ $(i) + $(um)*idy ] = ($(mb)[ $(i) + $(um)*idy ] - sum_i ) / $(mu)[ $(i+um*i) ] ;\n")         # bm = lm
        print(sb,"  $(mb)[ $(i) + $(um)*idx ] = ($(mb)[ $(i) + $(um)*idx ] - sum_i ) / $(mu)[ $(i+um*i) ] ;\n")         # bm = lm
        print(sb , "}\n")
        print(sb, "barrier(CLK_LOCAL_MEM_FENCE);\n")
    end
    print(sb,"// Solve U System DONE!\n")
    return String( take!( sb ) )
end



function ocl_print_matrix( name::String , sm::Int , sn::Int)
    return "if(idx==0&&idy==0){
              printf(\"\\n$(name) : $(sm)x$(sn) \\n\");
              for(int cm=0;cm<$(sm);cm++){ for( int cn=0;cn<$(sn);cn++){ printf( \" %6.3f \" , $(name)[cm+$(sm)*cn] ); }
                  printf(\"\\n\");
              }
              printf(\"\\n\");
             }
            "
end


function test_0()


    a = rand(Float32, 50_000)
    b = rand(Float32, 50_000)
    device, ctx, queue = OpenCL.cl.create_compute_context()

    a_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=a)
    b_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=b)
    c_buff = cl.Buffer(Float32, ctx, :w, length(a))

    p = cl.Program(ctx, source=test_source) |> cl.build!
    #cl.CL_device_info(1)

    k = cl.Kernel(p, "fwsim")
    queue(k, [2000], nothing, a_buff, b_buff, c_buff)
    r = cl.read(queue, c_buff)
end
