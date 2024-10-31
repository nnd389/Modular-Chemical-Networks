# %%
using Revise
using DrWatson
#@quickactivate "AstroChemNetwork"

    # %% Load packages 
using Catalyst
using DifferentialEquations
using Plots
using Sundials, LSODA
using OrdinaryDiffEq
using Symbolics
using DiffEqDevTools
using ODEInterface, ODEInterfaceDiffEq
using ModelingToolkit



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

    # %% Glover_network: Set up as a list of reactions. This is a network taken from a paper written by Glover
function run_Glover(d :: Dict)
    @unpack tspan, u0, params = d
    prob = ODEProblem(ssys, u0, tspan, params)
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
sol_NL99 = run_NL99(allvars)
#sol_Glover = run_Glover(allvars)

#sol_Glover = run_network(ssys, allvars)



    # %% Plot
#plot(sol_Benchmark, vars = (0,4), lw = 3)
#plot(sol_NL99_GC, vars = (0,16), lw = 3)
#plot(sol_NL99_odes, vars = (0,11), lw = 3, title = "HCO+ set up as a system of ODEs")
plot(sol_NL99, vars = (0,11), lw = 3, title = "HCO+ set up as a list of reactions")
#plot(sol_Glover, vars = (0,10), lw = 3, title = "HCO+ set up as a list of reactions")


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


