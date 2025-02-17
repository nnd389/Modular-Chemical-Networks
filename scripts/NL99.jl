
#using Revise
#using DrWatson
#@quickactivate "AstroChemNetwork"

    # %% Load packages 
using Catalyst
using DifferentialEquations
using Plots
using LSODA
using OrdinaryDiffEq
using Symbolics
using DiffEqDevTools
using ODEInterface, ODEInterfaceDiffEq
using ModelingToolkit


### Set The timespan, parameters, and initial conditions ###
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs

params = Dict(
    :T => 10, 
    :Av => 2, 
    :Go => 1.7, 
    :n_H => 611, 
    :shield => 1)

u0 = [0.5,    # 1:  H2   yep?
    9.059e-9, # 2:  H3+  yep
    2.0e-4,   # 3:  e    yep
    0.1,      # 4:  He  SEE lines 535 NL99
    7.866e-7, # 5:  He+  yep? should be 2.622e-5
    0.0,      # 6:  C    yep
    0.0,      # 7:  CHx  yep
    0.0004,   # 8:  O    yep
    0.0,      # 9:  OHx  yep
    0.0,      # 10: CO   yep
    0.0,      # 11: HCO+ yep
    0.0002,   # 12: C+   yep
    2.0e-7,   # 13: M+   yep
    2.0e-7]   # 14: M    yep


### Network admin things ###
#allvars = @strdict params tspan u0
@variables t 
@species H2(t) H3⁺(t) e(t) He(t) He⁺(t) C(t) CHx(t) O(t) OHx(t) CO(t) HCO⁺(t) C⁺(t) M⁺(t) M(t)
@parameters T Av Go n_H shield

### Define the network ###
reaction_equations = [
    # Cosmic-ray Ionization
    (@reaction 1.2e-17, H2 --> H3⁺ + e), # H2 --> H3⁺ + e + H
    (@reaction 6.8e-18, He --> He⁺ + e),

    # Ion-Molecule Reactions
    (@reaction n_H * 2e-9, H3⁺ + C --> CHx + H2),
    (@reaction n_H * 8e-10, H3⁺ + O --> OHx + H2),
    (@reaction n_H * 1.7e-9, H3⁺ + CO --> HCO⁺ + H2), 
    (@reaction n_H * 7e-15, He⁺ + H2 --> He), #FIX He⁺ + H2 --> He + H + H+
    (@reaction n_H * 1.6e-9, He⁺ + CO --> C⁺ + O + He),
    (@reaction n_H * 4e-16, C⁺ + H2 --> CHx), # C⁺ + H2 --> CHx + H
    (@reaction n_H * 1e-9, C⁺ + OHx --> HCO⁺),

    # Neutral-Neutral Reactions
    (@reaction n_H * 2e-10, O + CHx --> CO), # O + CHx --> CO + H
    (@reaction n_H * 5.8e-12 * (T)^(0.5), C + OHx --> CO), # C + OHx --> CO

    # Electron recombination
    (@reaction n_H * 9e-11 * (T)^(-0.64), He⁺ + e --> He),
    (@reaction n_H * 1.9e-6 * (T)^(-.54), H3⁺ + e --> H2), # H3⁺ + e --> H2 + H
    (@reaction n_H * 1.4e-10 * (T)^(-0.61), C⁺ + e --> C),
    (@reaction n_H * 3.3e-5 * (T)^(-1.0), HCO⁺ + e --> CO), # HCO⁺ + e --> CO + H
    (@reaction n_H * 3.8e-10 * (T)^(-0.65), M⁺ + e --> M),

    # Charge-Transfer Reactions
    (@reaction n_H * 2e-9, H3⁺ + M --> M⁺ + e + H2), 

    # Photoreactions
    (@reaction 3e-10 * Go * exp(-3 * Av), C --> C⁺ + e), 
    (@reaction 1e-9 * Go * exp(-1.5 * Av), CHx --> C), # CHx --> C + H
    (@reaction 10e-10 * (shield) * Go * exp(-3 *  Av), CO --> C + O), # FIX 
    (@reaction 5e-10 * Go * exp(-1.7 * Av), OHx --> O), # OHx --> O + H
    (@reaction 2e-10 * Go * exp(-1.9 * Av), M --> M⁺ + e),
    (@reaction 1.5e-10 * Go * exp(-2.5 * Av), HCO⁺ --> CO) # HCO⁺ --> CO + H
]


### Turn the Network into a system of ODEs and Timing ###
print("\n\nNelson Catalyst:")
@named system = ReactionSystem(reaction_equations, t)
#@named sys = ODESystem(reaction_equations, t) # this doesn't work but I have hope for it one day, see https://github.com/SciML/MethodOfLines.jl/issues/117

print("\nTime to convert:")
@time convert(ODESystem, complete(system))
sys = convert(ODESystem, complete(system))

#print("Time to simplify:")
#@time complete(sys)
#simplified_sys = complete(sys)
print("Time to complete:")
@time complete(sys)
completed_sys = complete(sys)


print("Time to create the completed problem:")
@time ODEProblem(completed_sys, u0, tspan, params)
prob_complete = ODEProblem(completed_sys, u0, tspan, params)

print("Time to solve the completed system with Rodas4(): ")
@time solve(prob_complete, lsoda(), saveat = 1e10);
sol_complete = solve(prob_complete, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e9)

#plot(sol, idxs = (0,6), lw = 3, lc = "blue")
#plot!(sol, idxs = (0,12), lw = 3, lc = "orange", title = "Nelson: Abundance of C and C+")


### Ensemble Problem ###
print("\nEnsemble timing: ")
prob = ODEProblem(completed_sys, u0, tspan, params)
#prob = ODEProblem(simplified_sys, u0, tspan, params)

function prob_func(prob, i, repeat)
    remake(prob, u0 = rand() * prob.u0)
end

print("\nTime to make (all the remakes of) the Ensemble Problem:")
@time EnsembleProblem(prob, prob_func = prob_func)
ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)

print("Time to solve the Ensemble Problem using EnsembleSerial")
@time solve(ensemble_prob,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSerial(); trajectories = 100)
#sim = solve(ensemble_prob,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSplitThreads(); trajectories = 100)

print("Time to solve the Ensemble Problem using EnsembleSplitThreads")
@time solve(ensemble_prob,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSplitThreads(); trajectories = 100)
#sim = solve(ensemble_prob,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSplitThreads(); trajectories = 100)

print("Time to solve the Ensemble Problem using EnsembleThreads")
@time solve(ensemble_prob,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleThreads(); trajectories = 100)
#sim = solve(ensemble_prob,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleThreads(); trajectories = 100)

print("Time to solve the Ensemble Problem using EnsembleDistributed")
@time solve(ensemble_prob,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleDistributed(), trajectories = 100)
#sim = solve(ensemble_prob,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleDistributed(), trajectories = 100)

#plot(sim, idxs = (0,6), linealpha = 1, lw = 3, title = "Nelson Catalyst: Ensemble problem for C and C+")
#plot!(sim, idxs = (0,12), linealpha = 0.4, lw = 3)


### We're gonna for loop and GPU this baddie###
# CAUTION! The second for loop always runs faster than the first regardless for some reason, I think it has to do with setting up the @time macro?
print("\nFor loop timing: ")
@time begin
for i in 1:100
    u0_rand = rand(14)
    prob_complete_for = ODEProblem(completed_sys, u0_rand, tspan, params)
    sol_complete_for = solve(prob_complete_for, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e9)
    #plot!(sol_for, idxs = (0,10), lw = 3, lc = "blue")
    #plot!(sol_for, idxs = (0,9), lw = 3, lc = "orange", title = "Glover with Glover rates: C and C+")
end
end

#=
### Plotting ### (ordering is off)

# C and C+
plot(sol, idxs = (0,6), lw = 3, lc = "blue")
plot!(sol, idxs = (0,12), lw = 3, lc = "orange", title = "Nelson: Abundance of C and C+")
#savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/C_and_Cp_Nelson.png")

# CO
plot(sol, idxs = (0,10), lw = 3, lc = "green", title = "Nelson: Abundance of CO")
#savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/CO_Nelson.png")

# O 
plot(sol, idxs = (0,8), lw = 3, lc = "blue", title = "Nelson: Abundance of O")
#savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/O_Nelson.png")

# He+
plot(sol, idxs = (0,5), lw = 3, lc = "light pink", title = "Nelson: Abundance of He+")
#savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/Hep_Nelson.png")

# CH and CH2 = CHx
plot(sol, idxs = (0,7), lw = 3, lc = "blue", title = "Nelson: CHx")
#savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/CH_and_CH2_Nelson.png")

# OH, OH+, H2O, H2O+, and O2 = OHx
plot(sol, idxs = (0,9), lw = 3, lc = "green", title = "Nelson: OHx")
#savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/OH_OHp_H2O_H2Op_O2_Nelson.png")

# H3+
plot(sol, idxs = (0,2), lw = 3, lc = "orange", title = "Nelson")
#savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/H3p_Nelson.png")

# HCO+ 
plot(sol, idxs = (0,11), lw = 3, lc = "orange", title = "Nelson")
#savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/HCOp_Nelson.png")

# M 
plot(sol, idxs = (0,14), lw = 3, lc = "light blue", title = "Nelson")
#savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/Mg_Fe_Na_Nelson.png")

=#





