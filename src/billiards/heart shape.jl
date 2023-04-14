function heart_eq(t) where T<:Real
    return SVector(16*sin(t)^3, 13 * cos(t) - 5*cos(2t) - 2*cos(3t) - cos(4*t))
end