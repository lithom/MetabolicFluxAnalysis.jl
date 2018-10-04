
struct rfd_parametrization
    p_raw::Polyhedron   # original polytope
    p_fd::Polyhedron    # full dimensional polytope
    p_opt::Polyhedron  # rounded, full dimensional polytope
    at_opt2raw::AffineMap
    at_opt2fd::AffineMap
    at_fd2raw::AffineMap

    cheby_p_fd::Vector{Float64}
    cheby_p_opt::Vector{Float64}
end

function eval_opt2raw( rfd_parametrization , x_opt::Array{Float64,2} )
    return rfd_parametrization.at_opt2raw(x_opt)
end

function eval_opt2raw( rfd_parametrization , x_opt::Array{Float64,1} )
    return rfd_parametrization.at_opt2raw(x_opt)
end


function get_rounded_fulldimensional_parametrization( fp::Polyhedron ; remove_red_constr=0)


    ( fp_fd , at_fd ) = MetabolicFluxAnalysis.find_minimal_affine_subspace(fp)
    if(remove_red_constr==1)
        fp_fd_a = MetabolicFluxAnalysis.remove_redundant_constraints( fp_fd )
        fp_fd = fp_fd_a
    end
    ( fp_rounded , at_rounded ) = round_polytope( fp_fd )
    at_opt2fd   = at_rounded
    at_opt2flux = compose(at_fd , at_rounded)

    cheby_p_opt  = zeros(size(fp_rounded.G,2))
    cheby_p_fd   = at_opt2fd(cheby_p_opt)

    return rfd_parametrization( copy(fp) , fp_fd , fp_rounded , at_opt2flux , at_opt2fd , at_fd , cheby_p_fd , cheby_p_opt  )
end
