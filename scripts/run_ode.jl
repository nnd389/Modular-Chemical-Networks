# %%
print("just checking 1\n")
using Revise
using DrWatson
#@quickactivate "AstroChemNetwork"

    # %% Load packages 
    print("just checking 2\n")
using Catalyst
using DifferentialEquations
using Plots
using Sundials, LSODA
using OrdinaryDiffEq
using Symbolics
using DiffEqDevTools
using ODEInterface, ODEInterfaceDiffEq
using ModelingToolkit
print("just checking 3\n")



    # %% File admin things
includet("Benchmark_network.jl")
includet("NL99_GC_network.jl")
includet("NL99_network_ODEs.jl")
includet("NL99_network.jl")
includet("Glover_network.jl")



    # %% Benchmark_network: This is a network taken from an astrochemistry benchmarking code
function run_Benchmark(d :: Dict)
    @unpack tspan, u0, params = d
    prob = ODEProblem(ssys, u0, tspan, params)
    sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    return sol
end

    # %% NL99_GC_network: This is a network taken from a paper written by Nelson and Langer in 1999 with some added hydrogen chemistry as combined in ____
# FIX: this network isn't complete, I just have to fill in some rates and double check
function run_NL99_GC(d :: Dict)
    @unpack tspan, u0, params = d
    prob = ODEProblem(ssys, u0, tspan, params)
    sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    return sol
end

    # %% NL99_network_ODEs: Set up as a system of ODEs instead of list of reactions, see NL99_network for set up as a list of reactions. This is a network taken from a paper written by Nelson and Langer in 1999
function run_NL99_odes(d :: Dict)
    @unpack tspan, u0, params = d
    prob = ODEProblem(NL99_network_odes, u0, tspan, params)
    sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    return sol
end

    # %% NL99_network: Set up as a list of reactions. This is a network taken from a paper written by Nelson and Langer in 1999 with some added hydrogen chemistry as combined in ____
function run_NL99(d :: Dict)
    @unpack tspan, u0, params = d
    prob = ODEProblem(ssys, u0, tspan, params)
    sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    return sol
end







k1 = 10^(-17.845 + 0.762 * log(T) + 0.1523 * (log(T))^2 - 0.03274 * (log(T))^3)
k2 = 0.0000000015
k6 = 0.00000001
k12 = 0.0000000000001269 * ((315614/T)^(1.503)) * (1 + (604625/T)^(0.47))^(-1.923) 
k14 = 0.0000000025634 * (T)^(1.78186) 
k15 = 0.0000000069 * T^(-0.35) 
k17 = 0.0000000001 * T^(-0.5) * (12.72 - 1.615 * log(T) - 0.3162 * (log(T))^2 + 0.0493 * (log(T))^3) 
k19 = 0.00000000126 * T^(-0.75) * exp(-127500/T) 
k20 = 0.00000000000467 * (T/300)^(-0.6) 
k21 = 0.00000000013 * T^(-0.64)
k29 = 0.0000000000000000858 * (T)^(0.757) 
k38 = 0.000000000066 
k42 = 0.00000000005 * (T/300)^(0.5)
k47 = 0.000000000035
k52 = 0.000000000047 * (T/300)^(-0.34)
k146 = 2.1e-19 
k149 = 2.5e-18
k154 = 1.32e-32 * (T/300)^-0.38 
k157 = 5.99e-33 * (T/5000)^(-1.6) 
k158 = 6.16e-29 * (T/300)^(-3.08) 
# FIX K12, K14, K17
# CHECK K22 and K23 and K30

R188 = 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2)
R189 = 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2)
R190 = 5.0e-10 * exp(-2.55*Av + 0.0165*(Av)^2)
R191 = 1.5e-10 * exp(-2.55*Av + 0.0165*(Av)^2)
R192 = 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2)
R193 = 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2)
R194 = 7.5e-12 * exp(-2.55*Av + 0.0165*(Av)^2)
R195 = 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2)

    # %% Glover_network: Set up as a list of reactions. This is a network taken from a paper written by Glover
function run_Glover(d :: Dict)
    
    print("just checking 4\n")
    @unpack tspan, u0 = d
    print("just checking 5\n")
    prob = ODEProblem(ssys, u0, tspan, params)
    print("just checking 6\n")
    sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    return sol
end




function run_network(ssys, d :: Dict)
    @unpack tspan, u0, params = d
    prob = ODEProblem(ssys, u0, tspan, params)
    sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    return sol
end







    # %% Run
#sol_Benchmark = run_Benchmark(allvars)
#sol_NL99_GC = run_NL99_GC(allvars)
#sol_NL99_odes = run_NL99_odes(allvars)
#sol_NL99 = run_NL99(allvars)
sol_Glover = run_Glover(allvars)

#sol_Glover = run_network(ssys, allvars)



    # %% Plot
#plot(sol_Benchmark, vars = (0,4), lw = 3)
#plot(sol_NL99_GC, vars = (0,16), lw = 3)
#plot(sol_NL99_odes, vars = (0,11), lw = 3, title = "HCO+ set up as a system of ODEs")
#plot(sol_NL99, vars = (0,11), lw = 3, title = "HCO+ set up as a list of reactions")
plot(sol_Glover, vars = (0,5), lw = 3, title = "HCO+ from Glover network")


#plot(sol_NL99, vars = (0,6), lw = 3, title = "C and C+ set up as a list of reactions")
#plot!(sol_NL99, vars = (0,12), lw = 3)
#plot(sol_NL99_odes, vars = (0,11), lw = 3, title = "HCO+ set up as a system of ODEs")

#plot(p1,p2)




    # THINGS TO FIX 
# JuliaGPU!!
# FIX GLOVER, check reaction rates, make conditional rates?, fix initial conditions, fix commented reactions, add n_H to two body reactions
# IMPORTANT: it looks like you have to run the function file before running this one 
# whenever you change networks the parameters are messing with each other
# Combine all functions into one
# verify results! check with DESPOTIC for all networks
# initial conditions for benchmark, NL99_GC are wrong
# Fix rates for NL_99_gc
# Might have to add n_H factor to the reaction rates of NL99_GC or any network that is defined as just the reactions not the ODES




    # %% Save results
#wsave(datadir("sims", savename(params)))


