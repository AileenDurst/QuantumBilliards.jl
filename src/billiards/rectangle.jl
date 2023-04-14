#include("../abstracttypes.jl")
#include("../basis/fourierbessel.jl")
#include("geometry.jl")
using StaticArrays
#this is old code and must be reworked

function make_quarter_rectangle(chi; h=one(chi),x0=zero(chi),y0=zero(chi),rot_angle=zero(chi))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    l = chi * h
    corners = [SVector(l/2, 0), SVector(l/2, h/2), SVector(zero(chi), h/2), SVector(zero(chi), zero(chi))]

    line1 = LineSegment(corners[1],corners[2];origin=origin,rot_angle=rot_angle)
    line2 = LineSegment(corners[2],corners[3];origin=origin,rot_angle=rot_angle)
    line3 = VirtualLineSegment(corners[3],corners[4];origin=origin,rot_angle=rot_angle)
    line4 = VirtualLineSegment(corners[4],corners[1];origin=origin,rot_angle=rot_angle)
    boundary = [line1, line2, line3, line4]
    return boundary, corners
end

function make_full_rectangle(chi;h=one(chi),x0=zero(chi),y0=zero(chi),rot_angle=zero(chi))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    l = chi * h
    corners = [SVector(l/2, -h/2), SVector(l/2, h/2), SVector(-l/2, h/2), SVector(l/2, -h/2)]

    line1 = LineSegment(corners[1],corners[2];origin=origin,rot_angle=rot_angle)
    line2 = LineSegment(corners[2],corners[3];origin=origin,rot_angle=rot_angle)
    line3 = LineSegment(corners[3],corners[4];origin=origin,rot_angle=rot_angle)
    line4 = LineSegment(corners[4],corners[1];origin=origin,rot_angle=rot_angle)
    boundary = [line1, line2, line3, line4]
    return boundary, corners
end



struct Rectangle{T}  <: AbsBilliard
    fundamental_boundary::Vector
    full_boundary::Vector
    length::T
    area::T
    chi::T
    h::T
    corners::Vector{SVector{2,T}}
end

function Rectangle(chi::Float64; curve_types = [:Virtual,:Real, :Real, :Virtual] , x0=0.0, y0=0.0, h=1.0)
    
    #println("α=$alpha, β=$beta, γ=$gamma")
    l = chi * h
    full_boundary, corners  = make_full_rectangle(chi;h=h,x0=x0,y0=y0)
    fundamental_boundary, _ = make_quarter_rectangle(chi;h=h,x0=x0,y0=y0)
    length = 2*l +2*h 
    area = h*l
    return Rectangle(fundamental_boundary,full_boundary,length,area,chi,h,corners)
end

function make_rectangle_and_basis(chi; curve_types = [:Virtual,:Real, :Real, :Virtual] , x0=0.0, y0=0.0, h=1.0)
    billiard = Rectangle(chi; curve_types = curve_types , x0=x0, y0=y0, h=h)
    basis = CornerAdaptedFourierBessel(1, pi/2.0, 0.0, x0, y0) 
    return billiard, basis 
end

