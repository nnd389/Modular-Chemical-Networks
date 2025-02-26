
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
#using XLSX
#using CSV, DataFrames


### Set The timespan, parameters, and initial conditions ###
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs

T = 10
Av = 2
Go = 1.7
n_H = 611
shield = 1

Grassi_T = 50
Grassi_ntot = 1e4
Grassi_zeta = 1e-16

params = Dict(
    :T => T, 
    :Av => Av, 
    :Go => Go, 
    :n_H => n_H, 
    :shield => shield)

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
@time convert(ODESystem, complete(system))
@time convert(ODESystem, complete(system))
sys = convert(ODESystem, complete(system))

#print("Time to simplify:")
#@time complete(sys)
#simplified_sys = complete(sys)
print("Time to complete:")
@time complete(sys)
@time complete(sys)
@time complete(sys)
completed_sys = complete(sys)


print("Time to create the completed problem:")
@time ODEProblem(completed_sys, u0, tspan, params)
@time ODEProblem(completed_sys, u0, tspan, params)
@time ODEProblem(completed_sys, u0, tspan, params)
prob_complete = ODEProblem(completed_sys, u0, tspan, params)

print("Time to solve the completed system with Rodas4(): ")
@time solve(prob_complete, lsoda(), saveat = 1e10);
@time solve(prob_complete, lsoda(), saveat = 1e10);
@time solve(prob_complete, lsoda(), saveat = 1e10);
sol_complete = solve(prob_complete, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e9)

plot(sol_complete, idxs = (0,6), lw = 3, lc = "blue")
plot!(sol_complete, idxs = (0,12), lw = 3, lc = "orange", title = "Nelson: Abundance of C and C+")


### Ensemble Problem ###
#=
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
=#




### We're gonna for loop and GPU this baddie###
num_runs = 3
plot()

for i in 1:num_runs
    u0_rand = rand().* u0
    ### Timespan
    spy = 3600 * 24 * 365
    years = 30000
    tstop = years * spy
    timestep = floor(tstop/100)
    tspan = (0, tstop)
    #timestep = 9.4608e9 # will spit out 101 steps for 30,000 years

    ### Solve
    println("\nDataset ", i)
    println("parameters use is: ", params)
    println("u0 used is: ", u0_rand)
    println("tspan used is: ", tspan)
    prob = ODEProblem(completed_sys, u0_rand, tspan, params)
    sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=timestep)
    plot!(sol, idxs = (0,6), lw = 1, lc = "orange")
    plot!(sol, idxs = (0,12), lw = 1, lc = "blue")
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





#=
### Creating Datasets ###
### Grassi's Parameters and Random Initialization
function Grassi_params()
    params_Grassi = Dict(:T => 10, :Av => 2, :Go => 1.7, :n_H => 611, :shield => 1)
    params_Grassi[:T] = 50
    params_Grassi[:n_H] = 10^4
    return params_Grassi
    # FIX to add zeta = cr_ion_rate
end

function Grassi_random_vector()
    vec = ones(14) 
    #u0 = [H2, H3+, e,  He, He+, C, CHx, O,  OHx, CO, HCO+, C+, M+, M]
    # I am following T. Grassi's CODE as close as possible, 
    # But it doesn't make sense to me with what the paper says
    ngas = 1e4
    rmax = -4
    rmin = -6

    vec[1] = ngas         # H2
    vec[2] = 1e-20 * ngas # H3+
    vec[3] = 4*(1e-20 * ngas) + ngas * 10^(rand() * (rmax-rmin) + rmin) # e = H3+ & He+ & HCO+ & C+ & M+ 
    vec[4] = 1e-20 * ngas # He
    vec[5] = 1e-20 * ngas # He+
    vec[6] = ngas * 10^(rand() * (rmax-rmin) + rmin) # C 
    vec[7] = 1e-20 * ngas # CHx
    vec[8] = 4*(ngas * 10^(rand() * (rmax-rmin) + rmin)) # O
    vec[9] = 1e-20 * ngas # OHx
    vec[10] = 1e-20 * ngas # CO
    vec[11] = 1e-20 * ngas # HCO+
    vec[12] = ngas * 10^(rand() * (rmax-rmin) + rmin)# C+
    vec[13] = 1e-20 * ngas # M+
    vec[14] = 1e-20 * ngas # M
    return vec
end



# CAUTION! The second for loop always runs faster than the first regardless for some reason, I think it has to do with setting up the @time macro?
num_runs = 3
plot()

for i in 1:num_runs
    ### Hydro Snapshot params + Nina initial conditions
    # Get the parameters from the hydro_snapshot .csv file
    # params = Dict(:T => 10, :Av => 2, :Go => 1.7, :n_H => 611, :shield => 1)
    #df = CSV.read("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/data/hydro_snapshot.csv", DataFrame)
    #snap_T = df[i, 2] 
    #println("Tempurature from snapshot is: ", snap_T) # (row 2, column 1)
    #params[:T] = snap_T
    #u0_rand = u0
    #u0_rand = rand().* u0

    ### Grassi params + Grassi initial conditons
    #local u0_rand = Grassi_random_vector()
    #local params_rand = Grassi_params()
    local u0_rand = Grassi_random_vector()
    local params_rand = params

    ### Timespan
    spy = 3600 * 24 * 365
    global years = 10
    tstop = years * spy
    timestep = floor(tstop/100)
    local tspan_rand = (0, tstop)
    #timestep = 9.4608e9 # will spit out 101 steps for 30,000 years

    ### Solve
    println("\nDataset ", i)
    println("parameters use is: ", params_rand)
    println("u0 used is: ", u0_rand)
    println("tspan used is: ", tspan_rand)
    local prob = ODEProblem(completed_sys, u0_rand, tspan_rand, params_rand)
    sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=timestep)
    plot!(sol, idxs = (0,6), lw = 1, lc = "orange")
    plot!(sol, idxs = (0,12), lw = 1, lc = "blue")

    # Create solution dataset
    data = permutedims(hcat(sol.u...))
    headers = vcat("Time (s)", "H2", "H3+", "e", "He", "He+", "C", "CHx", "O", "OHx", "CO", "HCO+", "C+", "M+", "M")
    save_dir = "/Users/kneenaugh/Desktop/Git/AstroChemNetwork/data"
    filename = joinpath(save_dir, "delete_run$(i).xlsx")

    XLSX.openxlsx(filename, mode="w") do xf
        sheet = xf[1]  # First sheet
        sheet["B2"] = data  # Write data 
        sheet["A1"] = headers  # Write headers
        sheet["A2"] = string(0) # Write time labels
        for i in 1:100
            sheet_str = "A" * string(i+2)
            t_label = string(timestep * i-1)
            sheet[sheet_str] = t_label
        end
        params_sheet = XLSX.addsheet!(xf, "Parameters")
        param_labels = vcat("T", "Av", "Go", "n_H", "shield")
        params_sheet["A1"] = param_labels
        params_sheet["A2"] = string(params[:T])
        params_sheet["B2"] = string(params[:Av])
        params_sheet["C2"] = string(params[:Go])
        params_sheet["D2"] = string(params[:n_H])
        params_sheet["E2"] = string(params[:shield])
    end
end

xlabel!("Time")
ylabel!("Concentration of H")
title = string("C and C+, Random u0 over ", years, " years")
title!(title)
display(plot!())  # Show the plot

=#



