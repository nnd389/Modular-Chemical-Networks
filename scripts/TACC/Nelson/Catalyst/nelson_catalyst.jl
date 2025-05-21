print("\n\n     Nelson Catalyst GPU Parrallelization Test     ")
print("\nCheck 1: We are at the begining of the Julia script!")
using Catalyst
using DifferentialEquations
using Plots
using LSODA
using OrdinaryDiffEq
using Symbolics
using DiffEqDevTools
using ODEInterface, ODEInterfaceDiffEq
using DiffEqGPU 
using CUDA
#using ModelingToolkit

print("\nCheck 2: Finished initializing packages")

### Set The timespan, parameters, and initial conditions ###
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs

params = [10,  # T
    2,   # Av
    1.7, # Go
    611, # n_H
    1]   # shield

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

 print("\nCheck 3: Finished initializing timespan, initial conditions, and parameters")


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
print("\nCheck 4: Finished reading the reactions equations")


### Turn the Network into a system of ODEs and Timing ###
@named system = ReactionSystem(reaction_equations, t)
#@named sys = ODESystem(reaction_equations, t) # this doesn't work but I have hope for it one day, see https://github.com/SciML/MethodOfLines.jl/issues/117
sys = convert(ODESystem, Catalyst.complete(system))
completed_sys = Catalyst.complete(sys)




print("\n\nNelson Catalyst CPUs Parrallelization Test:")
prob = ODEProblem(completed_sys, u0, tspan, params)
sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e9)

print("\nTime to solve the completed Nelson system ONCE with CPU lsoda: ")
@time solve(prob, lsoda(), saveat = 1e10);
@time solve(prob, lsoda(), saveat = 1e10);
@time solve(prob, lsoda(), saveat = 1e10);
@time solve(prob, lsoda(), saveat = 1e10);


### We're gonna for loop and GPU this baddie###
# CAUTION! The second for loop always runs faster than the first regardless for some reason, I think it has to do with setting up the @time macro?
num_runs = 10
print("\nFor loop timing for ", num_runs, " runs: ")

@time begin
    for i in 1:num_runs
        u0_rand = rand(Float32,14) .* u0
        prob_for = ODEProblem(completed_sys, u0_rand, tspan, params)
        sol_for = solve(prob_for, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    end
end
@time begin
    for i in 1:num_runs
        u0_rand = rand(Float32,14) .* u0
        prob_for = ODEProblem(completed_sys, u0_rand, tspan, params)
        sol_for = solve(prob_for, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    end
end
@time begin
    for i in 1:num_runs
        u0_rand = rand(Float32,14) .* u0
        prob_for = ODEProblem(completed_sys, u0_rand, tspan, params)
        sol_for = solve(prob_for, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    end
end
@time begin
    for i in 1:num_runs
        u0_rand = rand(Float32,14) .* u0
        prob_for = ODEProblem(completed_sys, u0_rand, tspan, params)
        sol_for = solve(prob_for, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
    end
end




### Ensemble Problem ###
print("\nEnsemble GPU Problem")
print("\nCheck 1: Finished the for loops")
tspan_ens = (0.0f0, 946080000000.0f0) # ~30 thousand yrs
print("\nCheck 2: Finished initializing the timespan")
prob_ens = ODEProblem(completed_sys, u0, tspan_ens, params) # perhaps try CUDA.cu(params) and CUDA.cu(u0)?
print("\nCheck 3: Finished creating the Ensemble ODE problem")
prob_func_nelson = (prob_ens, i, repeat) -> remake(prob_ens, u0 = rand(Float32,14) .* u0)
print("\nCheck 4: Finished making all the remakes")
monteprob_nelson = EnsembleProblem(prob, prob_func = prob_func_nelson, safetycopy = false);
print("\nCheck 5: Finished creating the Ensemble problem")


sol_ensemble = solve(monteprob_nelson, Rodas4(), EnsembleThreads(), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=10000000000.0f0)
print("\n\nEnsemble Timing to solve ", num_runs, " random Nelson systems on GPUs with Rodas4():")
@time solve(monteprob_nelson, Rodas4(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Rodas4(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Rodas4(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Rodas4(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)

print("\nEnsemble Timing to solve ", num_runs, " random Nelson systems on GPUs with Rodas5P():")
@time solve(monteprob_nelson, Rodas5P(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Rodas5P(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Rodas5P(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Rodas5P(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)

print("\nEnsemble Timing to solve ", num_runs, " random Nelson systems on GPUs with Tsit5():")
@time solve(monteprob_nelson, Tsit5(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Tsit5(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Tsit5(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Tsit5(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)

print("\nEnsemble Timing to solve ", num_runs, " random Nelson systems on GPUs with Vern9():")
@time solve(monteprob_nelson, Vern9(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Vern9(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Vern9(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(monteprob_nelson, Vern9(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = num_runs, reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)

print("\nWe're officially on for Nelson ODEs GPUs!! Onwards and march!\n")





### Some extra Ensemble Timing that I didn't feel like throwing away ###
print("\n\nSome extra Ensemble Timing that I didn't feel like throwing away")
print("Time to solve the Ensemble Problem using EnsembleSerial")
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSerial(); trajectories = 100)
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSerial(); trajectories = 100)
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSerial(); trajectories = 100)

print("Time to solve the Ensemble Problem using EnsembleSplitThreads")
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSplitThreads(); trajectories = 100)
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSplitThreads(); trajectories = 100)
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleSplitThreads(); trajectories = 100)

print("Time to solve the Ensemble Problem using EnsembleThreads")
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleThreads(); trajectories = 100)
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleThreads(); trajectories = 100)
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleThreads(); trajectories = 100)

print("Time to solve the Ensemble Problem using EnsembleDistributed")
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleDistributed(), trajectories = 100)
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleDistributed(), trajectories = 100)
@time solve(monteprob_nelson,lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10, EnsembleDistributed(), trajectories = 100)











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
