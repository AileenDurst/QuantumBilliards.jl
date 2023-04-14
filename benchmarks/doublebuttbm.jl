include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
#using Revise 
#using GLMakie
using CairoMakie
using StaticArrays
#using Latexify

b = 0.9;

billiard = DoubleButt(b)

basis_CAFB = CornerAdaptedFourierBessel(1, pi/2.0, SVector(0.,0.),0.) 

f = Figure(resolution = (1000,1000));
plot_geometry_test!(f, billiard)
display(f)

f = Figure(resolution = (1000,500))
plot_basis_test!(f, basis_CAFB, billiard)
display(f)


d = 5.0
b = 5.0
sw_solver = DecompositionMethod(d,b)
acc_solverA = ScalingMethodA(d,b)
acc_solverB = ScalingMethodB(d,b)
acc_solver = acc_solverA

k0 = 500.00
dk = 0.1
acc_infoA = benchmark_solver(acc_solverA, basis_CAFB, billiard, k0, dk; plot_matrix=true);
acc_infoB = benchmark_solver(acc_solverB, basis_CAFB, billiard, k0, dk; plot_matrix=true);

sw_info = benchmark_solver(sw_solver, basis_CAFB, billiard, k0, dk; plot_matrix=true, log=false);


f = Figure(resolution = (1000,500))
plot_solver_test!(f,sw_solver,basis_CAFB,billiard,5.0,10.0,0.01)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solverA,basis_CAFB,billiard,100.0,101.0,0.25)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solverB,basis_CAFB,billiard,100.0,101.0,0.25)
display(f)

f = Figure(resolution = (1000,500))
plot_solver_test!(f,acc_solver,basis_CAFB,billiard,500.0,501.0,0.125, tol = 1e-3)
display(f)




k0 = 100.0
dk = 0.1
k, ten = solve_wavenumber(acc_solver,basis_CAFB, billiard,k0,dk)
ks, ten = solve_spectrum(acc_solver,basis_CAFB, billiard,k0,dk)
k, ten = solve_wavenumber(sw_solver,basis_CAFB, billiard,k0,dk)
state = compute_eigenstate(sw_solver, basis_CAFB, billiard, k)

states = compute_eigenstate_bundle(acc_solver, basis_CAFB, billiard, k0;dk =0.1, tol=0.1)
states.X
states.ks
states.tens

f = Figure(resolution = (800,2500));
plot_probability!(f,state,b=10.0, log=(true,-5),inside_only=false)
display(f)

f = Figure(resolution = (1500,1500));
plot_momentum_function!(f,states;log=false)
display(f)

f = Figure(resolution = (1500,1500));
plot_husimi_function!(f,states;log=false)
display(f)

f = Figure(resolution = (1500,1500));
plot_state_test!(f,state; b_u= 10.0)
display(f)
