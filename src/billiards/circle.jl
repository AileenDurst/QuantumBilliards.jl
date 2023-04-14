
#include("../abstracttypes.jl")
#include("../basis/fourierbessel/corneradapted.jl")
#include("geometry.jl")

function make_quarter_circle(;radius=1.,x0=0.,y0=0.,rot_angle=0.)
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    circle = CircleSegment(radius,pi/2, 0. , 0., 0.; origin=origin, rot_angle = rot_angle)
    boundary = [circle]
    return boundary
end

function make_full_circle(;radius=1.,x0=0.,y0=0,rot_angle=0.)
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    circle = CircleSegment(radius,2*pi, 0.,0., 0.; origin=origin, rot_angle = rot_angle)
    boundary = [circle]
    return boundary
end

struct CircleBilliard{T}  <: AbsBilliard where {T<:Real}
    fundamental_boundary::Vector
    full_boundary::Vector
    length::T
    area::T
    radius::T
end

function CircleBilliard(radius=1.0,x0=0.0,y0=0.0)
    full_boundary = make_full_circle(radius=radius,x0=x0,y0=y0)
    area = (pi*radius^2)
    fundamental_boundary = make_quarter_circle(radius=radius,x0=x0,y0=y0)
    length = sum([crv.length for crv in full_boundary])
    #PolygonOps.area(collect(zip(x,y)))
    return CircleBilliard(fundamental_boundary,full_boundary,length,area,radius)
end 


