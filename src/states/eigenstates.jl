#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/typeutils.jl")

struct Eigenstate{K,T} <: StationaryState
    k::K
    vec::Vector{T}
    dim::Int
    eps::T
    #basis type
end

function Eigenstate(k, vec)  
    eps = set_precision(vec[1])
    if eltype(vec) <: Real
        filtered_vec = eltype(vec).([abs(v)>eps ? v : zero(vec[1]) for v in vec])
    else 
        filtered_vec = vec
    end
    return Eigenstate(k, filtered_vec, length(vec), eps)
end

function compute_eigenstate(solver::AbsSolver, basis::AbsBasis, billiard::AbsBilliard,k;sampler=gauss_legendre_nodes)
    L = real_length(billiard)
    dim = round(Int, L*k*solver.dim_scaling_factor/(2*pi))
    basis_new = resize_basis(basis, dim)
    pts = evaluate_points(solver, billiard, sampler, k)
    ten, vec = solve_vect(solver,basis_new, pts, k)
    return Eigenstate(k,vec)
end