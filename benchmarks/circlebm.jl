include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
#using Revise 
#using GLMakie
using CairoMakie
using StaticArrays
#using Latexify




billiard = CircleBilliard()

even_x_even_y = [XReflection(1), YReflection(1) , XYReflection(1,1)]
basis_CAFB = CornerAdaptedFourierBessel(1, pi/2.0, SVector(0.,0.),0.) 


f = Figure(resolution = (1000,1000))

plot_geometry_test!(f, billiard)
display(f)



d = 1.5
b = [5.0]
sw_solver   = DecompositionMethod(d,b)
acc_solverA = ScalingMethodA(d,b)
acc_solverB = ScalingMethodB(d,b)
acc_solver  = acc_solverA



k0 = 500.00
dk = 0.1


sw_info   = benchmark_solver(sw_solver, basis_CAFB, billiard, k0, dk; plot_matrix=true, log=false);
acc_infoA_CAFB = benchmark_solver(acc_solverA, basis_CAFB, billiard, k0, dk; plot_matrix=true);


f = Figure(resolution = (1000,500));
plot_solver_test!(f,sw_solver,basis_CAFB,billiard,5.0,10.0,0.01)
display(f)

f = Figure(resolution = (1000,500));
plot_solver_test!(f,acc_solver,basis_CAFB,billiard,100.0,101.0,0.15, tol = 1e-3)
display(f)




k0 = 302.0
dk = 0.1

k_CAFB, ten = solve_wavenumber(acc_solver, basis_CAFB, billiard,k0,dk)
state_CAFB = compute_eigenstate(sw_solver, basis_CAFB, billiard, k_CAFB)



f = Figure(resolution = (1500,1500));
plot_wavefunction!(f,state_CAFB; b= 5.0, dens = 100.0, fundamental_domain=false)
display(f)

f = Figure(resolution = (1500,1500));
plot_state_test!(f,state_CAFB; b_u= 10.0)
display(f)