using MKL
include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
#using Revise 
using CairoMakie
#using Latexify


gamma = 2/3*pi #sqrt(2)/2 * pi
chi  = 2.0

billiard, basis = make_triangle_and_basis(gamma, chi)

d = 3.5
b = 5.0
sw_solver = DecompositionMethod(d,b)
acc_solver = ScalingMethod(d,b)


k0 = 2000.001
dk = 0.01
acc_info = benchmark_solver(acc_solver, basis, billiard, gauss_legendre_nodes, k0, dk; plot_matrix=true);
sw_info = benchmark_solver(sw_solver, basis, billiard, gauss_legendre_nodes, k0, dk; plot_matrix=true, log=true);

f = Figure(resolution = (1000,500))
ax = Axis(f[1,1])
plot_geometry_test!(ax, billiard)
display(f)

f = Figure(resolution = (1000,500))
plot_basis_test!(f, basis, billiard)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,sw_solver,basis,billiard,5.0,10.0,0.01)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solver,basis,billiard,100.0,101.0,0.1)
display(f)

k0 = 100.001
dk = 0.01
k, ten = solve_wavenumber(acc_solver,basis, billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis, billiard, k)

f = Figure(resolution = (1500,1500))
plot_state_test!(f,state,basis,billiard)
display(f)


b_range =collect(range(2.0,6.0,step=0.5))
f = Figure(resolution = (1000,500))
plot_benchmarks!(f, sw_solver, basis, billiard, gauss_legendre_nodes, k0, dk, 3.5, b_range)
display(f)





