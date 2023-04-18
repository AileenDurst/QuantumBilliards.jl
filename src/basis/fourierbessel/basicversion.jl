using Bessels
using CoordinateTransformations, Rotations, StaticArrays







function BesselBasis_vector(r, phi, dim, k)
    part1 = Jv(0:dim+1, r*k)
   # part2 = zeros(dim)
    #for nu in 1:dim
    #    part2[nu] = iseven(nu) ? cos(nu*phi) : sin(nu*phi)
   # end
   # basis = part1[1:dim] .* part2
    
    deriv_r = zeros(dim)
    deriv_r[1] = - part1[2]
    for i in reverse(2:dim)
    
        deriv_r[i] = 0.5*(part1[i-1] - part1[i+1])
    end
    return part1, deriv_r
end


function BesselBasis_single(r, phi, nu, k)

    part2 = iseven(nu) ? cos(nu*phi) : sin(nu*phi)
    
    if nu == 0 
        part1 = Jv(0, r*k)
        deriv_r = -Jv(1, r*k)
    else
        part1 = Jv(nu-1:nu+1, r*k)
        deriv_r = 0.5*(part1[1] - part1[3])
    end
    basis = part1 .* part2
    return part1, deriv_r
end

#                   test of the basis 


rnage = [0,10,10]
Jv(nu, r) = Bessels.besselj(nu, r)  # nu is an integer

function Jvp(nu, r::T) where {T<:Real} # makes derivative of besselfunction
    let
    j_minus = Jv(nu-one(T),r)
    j_plus = Jv(nu+one(T),r)
    return 0.5*( j_minus - j_plus)
    end
end


rval = 1;
phival =0;
dimval = 3;
kval = 1;

basis1, deriv1 = BesselBasis_vector(rval,phival, dimval, kval);

basis2 = Jv(1, 1);
deriv2 = Jvp(1, 1);

basis3, deriv3 = BesselBasis_single(rval,phival, 1, kval);
print(deriv1, "      ", deriv2, "      ", deriv3)

print(basis1, "      ", basis2, "      ", basis3)

#                  test done


struct FourierBesselBasis{T} <: AbsBasis where  {T<:Real}
    cs::PolarCS{T}
    dim::Int64 #using concrete type
    #symmetries::Union{Vector{Sy},Nothing}
end


function FourierBesselBasis(dim::Int, origin, rot_angle)
    cs = PolarCS(origin, rot_angle)
    return FourierBesselBasis(cs, dim)
end



function resize_basis(basis::FourierBesselBasis, billiard::Bi, dim::Int, k) where {Bi<:AbsBilliard}
    if basis.dim == dim
        return basis
    else
        return FourierBessel(basis.cs, dim)
    end
end

function BesselBasis_single(r, phi, nu, k)
    part2 = iseven(nu) ? cos(nu*phi) : sin(nu*phi)
    if nu == 0 
        part1 = Jv(0, r*k)
        deriv_r = -Jv(1, r*k)
    else
        part1 = Jv(nu-1:nu+1, r*k)
        deriv_r = 0.5*(part1[1] - part1[3]) * part2
    end
    basis = part1 .* part2
    return basis, deriv_r
end

function BesselBasis_singleState(r, phi, index, k) # indeces checked
    nu = index -1
    part2 = iseven(nu) ? cos.(nu*phi) : sin.(nu*phi)
    part1 = Jv(nu, r*k)
    
    basis = part1 .* part2
    return basis
end

function BesselBasis_singleBessel(r, index, k) # indeces checked
    nu = index -1
    part1 = Jv(nu, r*k)
    
    basis = part1 
    return basis
end

function BesselBasis_singleDerivative(r, phi, idx, k)
    nu = idx -1
    part2 = iseven(nu) ? cos(nu*phi) : sin(nu*phi)
    if nu == 0 
        deriv_r = -Jv(1, r*k) * part2 *r
    else
        deriv_r = 0.5*(Jv(nu-1, r*k) - part1[nu+1]) * part2 *r
    end

    return deriv_r
end


@inline function basis_fun(basis::FourierBesselBasis{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map
        pt_pol = (cartesian_to_polar(pm(pt)) for pt in pts)

        #norm::T = one(T)/sqrt(basis.dim)
        return collect(BesselBasis_singleState(pt[1], pt[2], i, k) for pt in pt_pol)
    end
end

#=function make_Basis_vec(indices::AbstractArray, k::T, r, phi)
    first_nu = indices[1]-1
    last_nu = indices[end] -1

    part1 = Jv(first_nu:last_nu, r*k)
    part2 = zeros(dim)

    for nu in first_nu:last_nu
        part2[nu] = iseven(nu) ? cos(nu*phi) : sin(nu*phi)
    end
    basis = part1 .* part2
    return basis
end             
=#


@inline function basis_fun(basis::FourierBesselBasis{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, pts=pts

        pt_pol = (cartesian_to_polar(pm(pt)) for pt in pts)
 
        #norm::T = one(T)/sqrt(basis.dim)
        M =  length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            B[:,i] .= (BesselBasis_singleState(pt[1], pt[2], i, k) for pt in pt_pol)
        end
        return B 
    end
end


@inline function dk_fun(basis::FourierBesselBasis{T}, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    #translation of coordiante origin
    let pm = basis.cs.local_map, pts=pts
        pt_pol = [cartesian_to_polar(pm(pt)) for pt in pts]
        #norm::T = one(T)/sqrt(basis.dim)
        r   = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2) 
        dk = zeros(length(r))
        for b in 1:length(r)
            dk[b] = BesselBasis_singleDerivative(r[b], phi, i, k)
        end
       
        return dk
    end
end
    

@inline function dk_fun(basis::FourierBesselBasis{T}, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, nu=basis.nu, pts=pts
        pt_pol = [cartesian_to_polar(pm(pt)) for pt in pts]
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        M =  length(pts)
        N = length(indices)
        dB_dk = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            dB_dk[:,i] .= @. BesselBasis_singleDerivative(r, phi, i, k)
        end
        return dB_dk
    end
end



function gradient(basis::FourierBesselBasis, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    nu = i-1
    let pm = basis.cs.local_map, pts=pts
        pt_xy  = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy) #local cartesian coords
        #norm::T = one(T)/sqrt(basis.dim)
        r   = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X   = getindex.(pt_xy,1)
        Y   = getindex.(pt_xy,2)

        j1  = zeros(length(r))
        j2  = similar(j1)
        j   = similar(j1)

        for b in eachindex(r)
            j1[b] =  BesselBasis_singleBessel(r[b], i-1, k)
            j[b]  =  BesselBasis_singleBessel(r[b], i, k)
            j2[b] =  BesselBasis_singleBessel(r[b], i+1, k)
        end
        
        dj     = 0.5 .*(j1 .- j2) 
        angle  = iseven(nu) ? cos.(nu*phi) : sin.(nu*phi)
        dangle = iseven(nu) ? -sin.(nu*phi)*nu : cos.(nu*phi)*nu
     
        #println(size(s))
        dx = @. (dj*k*X/r*angle +  j* dangle *Y/r^2)
        dy = @. (dj*k*Y/r*angle -  j* dangle *X/r^2)

    return dx, dy
    end
end

function gradient(basis::FourierBesselBasis, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, pts=pts
        #local cartesian coords
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy)
        #norm::T = one(T)/sqrt(basis.dim)
        r   = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X   = getindex.(pt_xy,1)
        Y   = getindex.(pt_xy,2)
        M   = length(pts)
        N   = length(indices)
        j1  = zeros(length(r))
        j2  = similar(j1)
        j   = similar(j1)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)

        Threads.@threads for i in eachindex(indices)
            nu = i -1
            for b in 1:length(r)
                j1[b] =  BesselBasis_singleBessel(r[b], i-1, k)
                j[b]  =  BesselBasis_singleBessel(r[b], i, k)
                j2[b] =  BesselBasis_singleBessel(r[b], i+1, k)
            end
        #println(size(j))
            dj = 0.5 .*(j1 .- j2) 
            angle  = iseven(nu) ? cos.(nu*phi) : sin.(nu*phi)
            dangle = iseven(nu) ? -sin.(phi) : cos.(phi)

            #println(size(s))
            dB_dx[:,i] = @. (dj*k*X/r*angle - j*dangle * 1/Y * 1/(1+(X/Y)^2))
            dB_dy[:,i] = @. (dj*k*Y/r*angle - j*dangle * X/Y^2 * 1/(1+(X/Y)^2))
        end
        #println(size(s))
    return dB_dx, dB_dy
    end
end


function basis_and_gradient(basis::FourierBesselBasis, i::Int, k::T, pts::AbstractArray) where {T<:Real}
    nu = i-1
    let pm = basis.cs.local_map, pts=pts
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy) #local cartesian coords
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        j = BesselBasis_singleState(r, phi,nu , k)
        #println(size(j))
        dj = 0.5*(Jv(nu-1, r*k) - part1[nu+1]) 
        angle = iseven(nu) ? cos(nu*phi) : sin(nu*phi)
        dangle = iseven(nu) ? -sin(phi) : cos(phi)

        #println(size(s))
        bf = j * angle
        dx = @. (dj*k*x/r*angle -  j*dangle *1/Y* 1/(1+(X/Y)^2))
        dy = @. (dj*k*y/r*angle - j*dangle *X/Y^2* 1/(1+(X/Y)^2))
    return bf, dx, dy

    end
end


function basis_and_gradient(basis::FourierBesselBasis, indices::AbstractArray, k::T, pts::AbstractArray) where {T<:Real}
    let pm = basis.cs.local_map, pts=pts
        #local cartesian coords
        pt_xy = collect(pm(pt) for pt in pts)
        pt_pol = collect(cartesian_to_polar(pt) for pt in pt_xy)
        #norm::T = one(T)/sqrt(basis.dim)
        r = getindex.(pt_pol,1)
        phi = getindex.(pt_pol,2)
        X = getindex.(pt_xy,1)
        Y = getindex.(pt_xy,2)
        M = length(pts)
        N = length(indices)
        B = zeros(T,M,N)
        dB_dx = zeros(T,M,N)
        dB_dy = zeros(T,M,N)
        Threads.@threads for i in eachindex(indices)
            nu = i -1
            j = BesselBasis_singleState(r, phi,nu , k)
        #println(size(j))
            dj = 0.5*(Jv(nu-1, r*k) - part1[nu+1]) 
            angle = iseven(nu) ? cos(nu*phi) : sin(nu*phi)
            dangle = iseven(nu) ? -sin(phi) : cos(phi)

            #println(size(s))
            B[:,i] .= @. j*anlge
            dB_dx[:,i] = @. (dj*k*x/r*angle -  j*dangle * 1/Y * 1/(1+(X/Y)^2))
            dB_dy[:,i] = @. (dj*k*y/r*angle - j*dangle * X/Y^2 * 1/(1+(X/Y)^2))
        end
        #println(size(s))
    return B, dB_dx, dB_dy


    end
end



