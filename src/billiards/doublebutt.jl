function make_quarter_doublebutt(half_separation;radius=one(half_separation),x0=zero(half_separation),y0=zero(half_separation),rot_angle=zero(half_separation))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    type   = typeof(half_separation)
    theta  = pi - acos(half_separation/radius)
    circle = CircleSegment(radius,theta, zero(type), half_separation, zero(type); origin=origin, rot_angle = rot_angle)
    y_2    = radius*sin(theta)
    x_1    = radius+half_separation  
    corners = [SVector(x_1, zero(type)), SVector(zero(type), y_2), SVector(zero(type), zero(type))]
    
    line1 = VirtualLineSegment(corners[2],corners[3];origin=origin,rot_angle=rot_angle)
    line2 = VirtualLineSegment(corners[3],corners[1];origin=origin,rot_angle=rot_angle)
    boundary = [circle, line1, line2]
    return boundary, corners
end



function make_full_doublebutt(half_separation;radius=one(half_separation),x0=zero(half_separation),y0=zero(half_separation),rot_angle=zero(half_separation))
    #d(x, y, x0, y0, x1, y1) = @.((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)
    origin = SVector(x0,y0)
    type = typeof(half_separation)
    theta = acos(half_separation/radius)
    circle1 = CircleSegment(radius,2*theta, -pi/2, half_separation, zero(type); origin=origin, rot_angle = rot_angle)
    circle2 = CircleSegment(radius,2*theta, pi/2, - half_separation, zero(type); origin=origin, rot_angle = rot_angle)

    y_0 = radius*sin(theta)
    corners = [SVector(zero(type), -y_0), SVector(zero(type), y_0)]
    boundary = [circle1, circle2]
    return boundary, corners
end


struct DoubleButt{T}  <: AbsBilliard
    fundamental_boundary::Vector
    full_boundary::Vector
    length::T
    area::T
    half_separation::T
    theta::T
    corners::Vector{SVector{2,T}}
end

function DoubleButt(half_separation::Float64; radius=one(half_separation) , x0=0.0, y0=0.0, h=1.0)

    theta = acos(half_separation/radius)
    full_boundary, corners  = make_full_doublebutt(half_separation;radius=radius,x0=x0,y0=y0)
    fundamental_boundary, _ = make_quarter_doublebutt(half_separation;radius=radius,x0=x0,y0=y0)
    length = 2* (2*pi*radius- radius*theta)
    area   =  2*(pi*radius^2- radius^2/2 *(theta - sin(theta)))
    return Rectangle(fundamental_boundary,full_boundary,length,area, theta, half_separation,corners)
end