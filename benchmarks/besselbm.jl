include("../src/QuantumBillards2.jl")
using Revise
#using .QuantumBilliards
#using Revise 
#using GLMakie
using CairoMakie
using StaticArrays
#using Plots
b = 0.9;

basis = FourierBesselBasis(10, SVector(0.,0.), 0.) 
basis2 = CornerAdaptedFourierBessel(10,pi*1, SVector(0.,0.), 0.)
billiard = DoubleButt(b)
f = Figure(resolution = (1000,500));
plot_basis_test!(f, basis2, billiard; i = 1)
display(f)


x = collect(range(-2,2,300));
y = collect(range(-2,2,300));
r = [sqrt.(x[i]^2+y[j]^2) for i in 1:length(x) for j in 1:length(y)];
phi = [atan(x[i] ,y[j]) for i in 1:length(x) for j in 1:length(y)];
nu = 1

first = reshape([Bessels.besselj(nu,r[i]*10) .* sin(nu*phi[i]) for i in 1:length(r)],(length(x),length(y)));
f = Figure(resolution = (1000,500));
plot_heatmap_balaced!(f,x,y,first) 
display(f)


points = [SVector(x[i], y[j]) for i in 1:length(x) for j in 1:length(y)]

gradx, grady = gradient(basis, 2, 10, points)
gradxplot    = reshape(collect(gradx),(length(x),length(y)))
gradyplot    = reshape(collect(grady),(length(x),length(y)))

f = Figure(resolution = (1000,500));
plot_heatmap_balaced!(f,x,y,gradyplot) 
display(f)



