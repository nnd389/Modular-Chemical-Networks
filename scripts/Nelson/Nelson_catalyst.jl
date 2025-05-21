# Note from Nina: I implemented to Nelson network in two ways, the first as a system of ODEs, 
# and the second using Julia's Catalyst package which allows you to write the reactions and 
# Catalyst will convert them to ODE form for you. 
# The script below implements the Nelson network in Catalyst form. 
# If you wish the see the ODEs, please see the script: Nelson_network_fordocs.jl



### Packages ###
using Catalyst
using Plots
using LSODA

### Set The timespan, parameters, and initial conditions ###
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 1000000 * seconds_per_year) # ~30 thousand yrs

T = 300
Av = 2
Go = 1.7
n_H = 611
shield = 1

params = Dict(
    :T => T, 
    :Av => Av, 
    :Go => Go, 
    :n_H => n_H, 
    :shield => shield)

    #= These are old but I love them too much to get rid of them
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
=#
u0_fiducial = [0.5,# 1:  H2 
1e-8,     # 2:  H3+
0,      # 3:  e
0.1,      # 4:  He 
0.0,      # 5:  He+
5e-9,     # 6:  C      
0.0,      # 7:  CHx = CH and CH2
1e-8,     # 8:  O     
8.001e-5, # 9:  OHx = OH, H2O, O2 and their ions  Ilse has 8e-5 for H2O, and 1e-8 for O2, sum = .00008001
9.9e-5,   # 10: CO   
9e-9,     # 11: HCO+  
1e-9,     # 12: C+     
0.0,      # 13: M+ = Mg+, Fe+, Ca+, and Na+       Ilse has 1e-11 (Mg+), 1e-11 (Fe+), sum = 2e-11
0.0]      # 14: M = Mg, Fe, Ca, and Na

u0_neutral = [0.5, # 1:  H2
0.0,       # 2:  H3+
2e-8,      # 3:  e
0.1,       # 4:  He
0.0,       # 5:  He+
9.9015e-5, # 6:  C
0.0,       # 7:  CHx = CH and CH2
9.9019e-5,  # 8:  O 
8.001e-5,  # 9:  OHx = OH, H2O, O2 and their ions
0.0,       # 10: CO
0.0,       # 11: HCO+
0.0,       # 12: C+
0.0,       # 13: M+ = Mg+, Fe+, Ca+, and Na+      Ilse has 1e-11 (Mg+), 1e-11 (Fe+), sum = 2e-11
0]         # 14: M = Mg, Fe, Ca, and Na 



### Network admin things ###
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



### Turn the Network into a system of ODEs ###
print("\n\nNelson Catalyst:")
@named system = ReactionSystem(reaction_equations, t)
#@named sys = ODESystem(reaction_equations, t) # this doesn't work but I have hope for it one day, see https://github.com/SciML/MethodOfLines.jl/issues/117

sys_nelson = convert(ODESystem, complete(system))
comp_nelson = complete(sys_nelson)
nelson_colors = palette(:devon, 10)


prob_nelson_fiducial = ODEProblem(comp_nelson, u0_fiducial, tspan, params)
sol_nelson_fiducial = solve(prob_nelson_fiducial, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=9.4608e10)
using Plots
Plots.plot(sol_nelson_fiducial, vars = (0,12), lw = 3, lc = nelson_colors[1], label = "fiducial", title = "C+")

prob_nelson_neutral = ODEProblem(comp_nelson, u0_neutral, tspan, params)
sol_nelson_neutral = solve(prob_nelson_neutral, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=9.4608e10)
Plots.plot!(sol_nelson_neutral, vars = (0,12), lw = 3, lc = nelson_colors[4], label = "neutral")


### Exporting to excel for plotting ###
#=
using DataFrames
using XLSX

function solution_to_df(sol, species_names)
    t = sol.t
    u = reduce(hcat, sol.u)'  # rows = time steps, columns = species
    df = DataFrame(time = t)
    for (i, name) in enumerate(species_names)
        df[!, name] = u[:, i]
    end
    return df
end

species_names = ["H2", "H3+", "e", "He", "He+", "C", "CHx", "O", "OHx", "CO", "HCO+", "C+", "M+", "M"]

df_fiducial = solution_to_df(sol_nelson_fiducial, species_names)
df_neutral = solution_to_df(sol_nelson_neutral, species_names)

XLSX.openxlsx("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/Data/nelson_comparison.xlsx", mode="w") do xf
    # Create and write to "fiducial Initial" sheet
    sheet1 = XLSX.addsheet!(xf, "fiducial Initial")
    XLSX.writetable!(sheet1, collect(eachcol(df_fiducial)), names(df_fiducial))

    # Create and write to "Neutral Initial" sheet
    sheet2 = XLSX.addsheet!(xf, "Neutral Initial")
    XLSX.writetable!(sheet2, collect(eachcol(df_neutral)), names(df_neutral))
end
=#



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
