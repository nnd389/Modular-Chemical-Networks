using OrdinaryDiffEq
using DiffEqDevTools, Plots
using LSODA

#=
The ODE function defined below models the reduced carbon-oxygen 
chemistry network of Nelson & Langer (1999, ApJ, 524, 923).

This Julia ODE function was written by Nina De La Torre advised by Dr. Stella Offner.
The solution was compared with results derived by DESPOTIC, (Mark Krumholz, 2013) 
a code to Derive the Energetics and Spectra of Optically Thick Insterstellar Clouds. 
DESPOTIC has pre-defined networks, one of them coming from the Nelson & Langer paper, 
so the initial conditions and parameters were meant to mimic those from DESPOTIC.

Note: The composite hydrocarbon radical CHx represents both CH and CH2, 
the composite oxygen species OHx represents OH, H2O, O2 and their ions,
and M represents the low ionization potential metals Mg, Fe, Ca, and Na. 

    Parameter definitions:
    T = 10     --> Tempurature (Kelvin)
    Av = 2     --> V-Band Extinction
    G₀ = 1.7   --> Go; "a factor that determines the flux of FUV radiation relative to the standard interstellar value (G₀ = 1) as reported by Habing (1968)."
    n_H = 611  --> Hydrogen Number Density
    shield = 1 --> "CO self-shielding factor of van Dishoeck & Black (1988), taken from Bergin et al. (1995)"
=#

function Nelson!(du,u,p,t)
    T, Av, Go, n_H, shield = p

    # 1: H2
    du[1] = -1.2e-17 * u[1] + 
            n_H * (1.9e-6 * u[2] * u[3]) / (T^0.54) - 
            n_H * 4e-16 * u[1] * u[12] - 
            n_H * 7e-15 * u[1] * u[5] + 
            n_H * 1.7e-9 * u[10] * u[2] + 
            n_H * 2e-9 * u[2] * u[6] + 
            n_H * 2e-9 * u[2] * u[14] + 
            n_H * 8e-10 * u[2] * u[8] 

    # 2: H3+
    du[2] = 1.2e-17 * u[1] + 
            n_H * (-1.9e-6 * u[3] * u[2]) / (T^0.54) - 
            n_H * 1.7e-9 * u[10] * u[2] - 
            n_H * 2e-9 * u[2] * u[6] - 
            n_H * 2e-9 * u[2] * u[14] - 
            n_H * 8e-10 * u[2] * u[8]

    # 3: e
    du[3] = n_H * (-1.4e-10 * u[3] * u[12]) / (T^0.61) - 
            n_H * (3.8e-10 * u[13] * u[3]) / (T^0.65) - 
            n_H * (3.3e-5 * u[11] * u[3]) / T + 
            1.2e-17 * u[1] - 
            n_H * (1.9e-6 * u[3] * u[2]) / (T^0.54) + 
            6.8e-18 * u[4] - 
            n_H * (9e-11 * u[3] * u[5]) / (T^0.64) + 
            3e-10 * Go * exp(-3 * Av) * u[6] +
            n_H * 2e-9 * u[2] * u[13]
            + 2.0e-10 * Go * exp(-1.9 * Av) * u[14]

    # 4: He
    du[4] = n_H * (9e-11 * u[3] * u[5]) / (T^0.64) - 
            6.8e-18 * u[4] + 
            n_H * 7e-15 * u[1] * u[5] + 
            n_H * 1.6e-9 * u[10] * u[5]

    # 5: He+
    du[5] = 6.8e-18 * u[4] - 
            n_H * (9e-11 * u[3] * u[5]) / (T^0.64) - 
            n_H * 7e-15 * u[1] * u[5] - 
            n_H * 1.6e-9 * u[10] * u[5]

    # 6: C
    du[6] = n_H * (1.4e-10 * u[3] * u[12]) / (T^0.61) - 
            n_H * 2e-9 * u[2] * u[6] - 
            n_H * 5.8e-12 * (T^0.5) * u[9] * u[6] + 
            1e-9 * Go * exp(-1.5 * Av) * u[7] - 
            3e-10 * Go * exp(-3 * Av) * u[6] + 
            1e-10 * Go * exp(-3 * Av) * u[10] * shield

    # 7: CHx
    du[7] = n_H * (-2e-10) * u[7] * u[8] + 
            n_H * 4e-16 * u[1] * u[12] + 
            n_H * 2e-9 * u[2] * u[6] - 
            1e-9 * Go * u[7] * exp(-1.5 * Av)

    # 8: O
    du[8] = n_H * (-2e-10) * u[7] * u[8] + 
            n_H * 1.6e-9 * u[10] * u[5] - 
            n_H * 8e-10 * u[2] * u[8] + 
            5e-10 * Go * exp(-1.7 * Av) * u[9] + 
            1e-10 * Go * exp(-3 * Av) * u[10] * shield

    # 9: OHx
    du[9] = n_H * (-1e-9) * u[9] * u[12] + 
            n_H * 8e-10 * u[2] * u[8] - 
            n_H * 5.8e-12 * (T^0.5) * u[9] * u[6] - 
            5e-10 * Go * exp(-1.7 * Av) * u[9]

    # 10: CO
    du[10] = n_H * (3.3e-5 * u[11] * u[3]) / T + 
             n_H * 2e-10 * u[7] * u[8] - 
             n_H * 1.7e-9 * u[10] * u[2] - 
             n_H * 1.6e-9 * u[10] * u[5] + 
             n_H * 5.8e-12 * (T^0.5) * u[9] * u[6] - 
             1e-10 * Go * exp(-3 * Av) * u[10] + 
             1.5e-10 * Go * exp(-2.5 * Av) * u[11] * shield

    # 11: HCO+
    du[11] = n_H * (-3.3e-5 * u[11] * u[3]) / T + 
             n_H * 1e-9 * u[9] * u[12] + 
             n_H * 1.7e-9 * u[10] * u[2] - 
             1.5e-10 * Go * exp(-2.5 * Av) * u[11]

    # 12: C+
    du[12] = n_H * (-1.4e-10 * u[3] * u[12]) / (T^0.61) - 
             n_H * 4e-16 * u[1] * u[12] - 
             n_H * 1e-9 * u[9] * u[12] + 
             n_H * 1.6e-9 * u[10] * u[5] + 
             3e-10 * Go * exp(-3 * Av) * u[6]

    # 13: M+
    du[13] = n_H * (-3.8e-10 * u[13] * u[3]) / (T^0.65) + 
             n_H * 2e-9 * u[2] * u[14] +
             2.0e-10 * Go * exp(-1.9 * Av) * u[14]

    # 14: M
    du[14] = n_H * (3.8e-10 * u[13] * u[3]) / (T^0.65) - 
             n_H * 2e-9 * u[2] * u[14] -
             2.0e-10 * Go * exp(-1.9 * Av) * u[14]

end

# Set the Timespan, Parameters, and Initial Conditions
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 1000000 * seconds_per_year) # ~30 thousand yrs

params = [10,  # T 10
          2,   # Av 2
          1.7, # Go 1.7
          611, # n_H 611
          1]   # shield 1

u0_fiducial = [0.5,# 1:  H2 
1e-8,     # 2:  H3+
0.0,      # 3:  e
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
9.9095e-5, # 6:  C
0.0,       # 7:  CHx = CH and CH2
1.79044e-4,# 8:  O 
8.001e-5,  # 9:  OHx = OH, H2O, O2 and their ions
0.0,       # 10: CO
0.0,       # 11: HCO+
0.0,       # 12: C+
0.0,       # 13: M+ = Mg+, Fe+, Ca+, and Na+      Ilse has 1e-11 (Mg+), 1e-11 (Fe+), sum = 2e-11
0]         # 14: M = Mg, Fe, Ca, and Na 



prob_neutral = ODEProblem(Nelson!, u0_neutral, tspan, params)
sol_neutral = solve(prob_neutral, Rodas4(), saveat = 1e10)

prob_fiducial = ODEProblem(Nelson!, u0_fiducial, tspan, params)
sol_fiducial = solve(prob_fiducial, Rodas4(), saveat = 1e10)


# Plotting
# C+
plot(sol_fiducial, vars = (0,12), lw = 3, lc = "orange", label = "C+ with fiducial u0", title = "Nelson C+ Comparision")
plot!(sol_neutral, vars = (0,12), lw = 3, lc = "green", label = "C+ with Neutral u0", xlabel = "1e7 years")

#=
# O
plot(sol_fiducial, vars = (0,8), lw = 3, lc = "blue", label = "O with fiducial u0", title = "Nelson O Comparision")
plot!(sol_neutral, vars = (0,8), lw = 3, lc = "green", label = "O with Neutral u0", xlabel = "1e6 years")

# CO
plot(sol_fiducial, vars = (0,10), lw = 3, lc = "blue", label = "CO with fiducial u0", title = "Nelson CO Comparision")
plot!(sol_neutral, vars = (0,10), lw = 3, lc = "green", label = "CO with Neutral u0", xlabel = "1e7 years")

# He
plot(sol_fiducial, vars = (0,4), lw = 3, lc = "blue", label = "He with fiducial u0", title = "Nelson He Comparision")
plot!(sol_neutral, vars = (0,4), lw = 3, lc = "green", label = "He with Neutral u0", xlabel = "1e9 years")

# e
plot(sol_fiducial, vars = (0,3), lw = 3, lc = "blue", label = "e with fiducial u0", title = "Nelson e Comparision")
plot!(sol_neutral, vars = (0,3), lw = 3, lc = "green", label = "e with Neutral u0", xlabel = "2e6 years")

# H3+
plot(sol_fiducial, vars = (0,2), lw = 3, lc = "blue", label = "H3+ with fiducial u0", title = "Nelson H3+ Comparision")
plot!(sol_neutral, vars = (0,2), lw = 3, lc = "green", label = "H3+ with Neutral u0", xlabel = "1e6 years")

# HCO+
plot(sol_fiducial, vars = (0,11), lw = 3, lc = "blue", label = "HCO+ with fiducial u0", title = "Nelson HCO+ Comparision")
plot!(sol_neutral, vars = (0,11), lw = 3, lc = "green", label = "HCO+ with Neutral u0", xlabel = "1e5 years")

=#





#=
colors = palette(:acton, 5)
p1 = plot(sol1, vars = (0,11), lc=colors[1], legend = false, titlefontsize = 12, lw = 3, xlabel = "", title = "HCO+ solved using Rodas5 (correct solution)")
p2 = plot(sol2, vars = (0,11), lc=colors[2], legend = false, titlefontsize = 12, lw = 3, xlabel = "", title = "HCO+ solved using FBDF")
p3 = plot(sol3, vars = (0,11), lc=colors[3], legend = false, titlefontsize = 12, lw = 3, xlabel = "", title = "HCO+ solved using lsoda")
p4 = plot(sol4, vars = (0,11), lc=colors[4], legend = false, titlefontsize = 12, lw = 3, xlabel = "", title = "HCO+ solved using lsoda with saveat")
combined_plot = plot(p1, p2, p3, p4, layout=(4, 1), dpi = 600, pallete=:acton)


@time solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
@time solve(prob, Vern9(), abstol=1e-14, reltol=1e-14)
@time solve(prob, Euler(), dt = 1e6)
=#
#=



# Run Benchmark
setups = [
          Dict(:alg=>Rodas4P()),
          Dict(:alg=>Rodas4()),
          Dict(:alg=>Rodas5P()),

          Dict(:alg=>FBDF()),
          Dict(:alg=>QNDF()),
          Dict(:alg=>CVODE_BDF()),
          Dict(:alg=>lsoda()),
          #Dict(:alg=>AutoTsit5(Rosenbrock23())),
          #Dict(:alg=>Euler(), :dt=>1e6)
          ]

refsol = solve(prob, Euler(); dt = 1e6)

abstols = 1.0 ./ 10.0 .^ (7:13)
reltols = 1.0 ./ 10.0 .^ (4:10)
#wp = WorkPrecision(prob,Rodas4(),appxsol=refsol,abstols,reltols;name="Dormand-Prince 4/5")
#plot(wp)

wp = WorkPrecisionSet(prob,
                      abstols,
                      reltols,
                      setups;
                      appxsol=refsol,
                      save_everystep=false, 
                      print_names = true)

plot(wp; 
     palette=:PRGn_8, 
     xlabel = "Error", 
     ylabel = "Time to solve (s)",
     xguidefontsize=16, 
     yguidefontsize=16, 
     tickfontsize=11)


#palette=:acton10
=#


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

df_fiducial = solution_to_df(sol_fiducial, species_names)
df_neutral = solution_to_df(sol_neutral, species_names)

XLSX.openxlsx("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/Data/nelson_comparison.xlsx", mode="w") do xf
    # Create and write to "Fiducial Initial" sheet
    sheet1 = XLSX.addsheet!(xf, "Fiducial Initial")
    XLSX.writetable!(sheet1, collect(eachcol(df_fiducial)), names(df_fiducial))

    # Create and write to "Neutral Initial" sheet
    sheet2 = XLSX.addsheet!(xf, "Neutral Initial")
    XLSX.writetable!(sheet2, collect(eachcol(df_neutral)), names(df_neutral))
end
=#


