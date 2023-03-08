#include("../abstracttypes.jl")
#include("../utils/billiardutils.jl")
#include("../utils/typeutils.jl")

struct Eigenstate{K,T} <: StationaryState
    k::K
    k_basis::K
    vec::Vector{T}
    ten::T
    dim::Int64
    eps::T
    #basis type
end

function Eigenstate(k, vec, ten)  
    eps = set_precision(vec[1])
    if eltype(vec) <: Real
        filtered_vec = eltype(vec).([abs(v)>eps ? v : zero(vec[1]) for v in vec])
    else 
        filtered_vec = vec
    end
    return Eigenstate(k, k, filtered_vec,ten, length(vec), eps)
end

function Eigenstate(k, k_basis, vec, ten)  
    eps = set_precision(vec[1])
    if eltype(vec) <: Real
        filtered_vec = eltype(vec).([abs(v)>eps ? v : zero(vec[1]) for v in vec])
    else 
        filtered_vec = vec
    end
    return Eigenstate(k, k_basis, filtered_vec, ten, length(vec), eps)
end

function compute_eigenstate(solver::SweepSolver, basis::AbsBasis, billiard::AbsBilliard,k;sampler=gauss_legendre_nodes)
    L = real_length(billiard)
    dim = round(Int, L*k*solver.dim_scaling_factor/(2*pi))
    basis_new = resize_basis(basis, dim)
    pts = evaluate_points(solver, billiard, sampler, k)
    ten, vec = solve_vect(solver,basis_new, pts, k)
    return Eigenstate(k, vec, ten)
end

function compute_eigenstate(solver::AcceleratedSolver, basis::AbsBasis, billiard::AbsBilliard, k;sampler=gauss_legendre_nodes, dk = 0.1)
    L = real_length(billiard)
    dim = round(Int, L*k*solver.dim_scaling_factor/(2*pi))
    basis_new = resize_basis(basis, dim)
    pts = evaluate_points(solver, billiard, sampler, k)
    ks, tens, X = solve_vectors(solver,basis_new, pts, k, dk)
    idx = findmin(abs.(ks.-k))[2]
    k_state = ks[idx]
    ten = tens[idx]
    vec = X[:,idx]
    return Eigenstate(k_state, k, vec, ten)
end

struct EigenstateBundle{K,T} <: AbsState
    ks::Vector{K}
    k_basis::K
    X::Matrix{T}
    tens::Vector{T}
    dim::Int64
    eps::T
    #basis type
end

function EigenstateBundle(ks, k_basis, X, tens)  
    eps = set_precision(X[1,1])
    type = eltype(X)
    if  type <: Real
        filtered_array = type.([abs(x)>eps ? x : zero(type) for x in X])
    else 
        filtered_array = X
    end
    return EigenstateBundle(ks, k_basis, filtered_array, tens, length(X[:,1]), eps)
end

function compute_eigenstate_bundle(solver::AcceleratedSolver, basis::AbsBasis, billiard::AbsBilliard, k;sampler=gauss_legendre_nodes, dk = 0.1, tol=1e-5)
    L = real_length(billiard)
    dim = round(Int, L*k*solver.dim_scaling_factor/(2*pi))
    basis_new = resize_basis(basis, dim)
    pts = evaluate_points(solver, billiard, sampler, k)
    ks, tens, X = solve_vectors(solver,basis_new, pts, k, dk)
    idx = abs.(tens) .< tol
    ks = ks[idx]
    tens = tens[idx]
    X = X[:,idx]
    return EigenstateBundle(ks, k, X, tens)
end