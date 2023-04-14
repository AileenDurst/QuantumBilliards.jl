include("../src/QuantumBilliards.jl")
#using Revise
using .QuantumBilliards
#using Revise 
#using GLMakie
using CairoMakie
using StaticArrays
#using Latexify


chi = 1.5;

billiard = Rectangle(chi)

basis_CAFB = CornerAdaptedFourierBessel(1, pi/2.0, SVector(0.,0.),0.) 
even_x_even_y = [XReflection(1), YReflection(1) , XYReflection(1,1)]

basis_RPW =  RealPlaneWaves(10, even_x_even_y;angle_arc=pi/2.0)

f = Figure(resolution = (1000,1000))
plot_geometry_test!(f, billiard)
display(f)



d = 2.
b = [5.0]
sw_solver   = DecompositionMethod(d,b)
acc_solverA = ScalingMethodA(d,b)
acc_solverB = ScalingMethodB(d,b)
acc_solver  = acc_solverA



k0 = 500.00
dk = 0.3


sw_info   = benchmark_solver(sw_solver, basis_CAFB, billiard, k0, dk; plot_matrix=true, log=false);
acc_infoA = benchmark_solver(acc_solverA, basis_CAFB, billiard, k0, dk; plot_matrix=true);

sw_info_RPW    = benchmark_solver(sw_solver, basis_RPW , billiard, k0, dk; plot_matrix=true, log=false);
acc_infoA_RPW  = benchmark_solver(acc_solverA, basis_RPW , billiard, k0, dk; plot_matrix=true);

f = Figure(resolution = (1000,500));
plot_solver_test!(f,sw_solver,basis_CAFB,billiard,5.0,10.0,0.01)
display(f)

f = Figure(resolution = (1000,500));
plot_solver_test!(f,acc_solver,basis_CAFB,billiard,100.0,101.0,0.5, tol = 1e-3)
display(f)


f = Figure(resolution = (1000,500));
plot_solver_test!(f,sw_solver,basis_RPW,billiard,0.5,4.0,0.01)
display(f)

f = Figure(resolution = (1000,500));
plot_solver_test!(f,acc_solver,basis_RPW,billiard,100.0,101.0,0.5, tol = 1e-3)
display(f)



k0 = 302.0
dk = 0.3

k_CAFB, ten = solve_wavenumber(acc_solver, basis_CAFB, billiard,k0,dk)
state_CAFB = compute_eigenstate(sw_solver, basis_CAFB, billiard, k_CAFB)

f = Figure(resolution = (1500,1500));
plot_wavefunction!(f,state_CAFB; b= 5.0, dens = 100.0, fundamental_domain=false)
display(f)

f = Figure(resolution = (1500,1500));
plot_state_test!(f,state_CAFB; b_u= 10.0)
display(f)


k_RPW, ten = solve_wavenumber(acc_solver, basis_RPW, billiard,k0,dk)
k_RPW, ten = solve_wavenumber(sw_solver, basis_RPW, billiard,3.8,0.4)
state_RPW= compute_eigenstate(sw_solver, basis_RPW, billiard, k_RPW)



f = Figure(resolution = (1500,1500));
plot_wavefunction!(f,state_RPW; b= 5.0, dens = 100.0, fundamental_domain=false)
display(f)

f = Figure(resolution = (1500,1500));
plot_state_test!(f,state_RPW; b_u= 10.0)
display(f)