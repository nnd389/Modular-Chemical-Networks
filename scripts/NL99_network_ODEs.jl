#=
The ODE function defined below models the reduced carbon-oxygen 
chemistry network of Nelson & Langer (1999, ApJ, 524, 923).

This Julia function was written by Nina De La Torre (email: haninadlt@gmail.com) advised by Dr. Stella Offner.
and the solution was compared with results derived by DESPOTIC,  (Mark Krumholz, 2013) 
a code to Derive the Energetics and Spectra of Optically Thick Insterstellar Clouds. 
DESPOTIC has pre-defined networks, one of them coming from the Nelson & Langer paper, 
so the initial conditions and paramters were meant to mimic those from DESPOTIC.

        Parameter definitions:
        T: Tempurature (Kelvin)
        Av: V-band extinction
        G₀: Or Go here; "a factor that determines the flux of FUV radiation relative to the standard interstellar value (G₀ = 1) as reported by Habing (1968)."
        n_H: Number density
        shield: "CO self-shielding factor of van Dishoeck & Black (1988), taken from Bergin et al. (1995)"

=#
    
    # %% Set The timespan, paramters, and initial conditions
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs

u0 = [0.5,      # 1:  H2   yep?
      9.059e-9, # 2:  H3+  yep
      2.0e-4,   # 3:  e    yep
      0.1,      # 4:  He   SEE lines 535 NL99
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

params = (10,  # T
          2,   # Av
          1.7, # Go
          611, # n_H
          1)   # shield
allvars = @strdict params tspan u0


# %% Nelson Network ODE definition
#=========================================#

function NL99_network_odes(du,u,p,t)
    T, Av, Go, n_H, shield = p
# 1: H2
#= du[1] = -1.2e-17 * u[1] + 
        n_H * (1.9e-6 * u[2] * u[3]) / (T^0.54) - 
        n_H * 4e-16 * u[1] * u[12] - 
        n_H * 7e-15 * u[1] * u[5] + 
        n_H * 1.7e-9 * u[10] * u[2] + 
        n_H * 2e-9 * u[2] * u[6] + 
        n_H * 2e-9 * u[2] * u[14] + 
        n_H * 8e-10 * u[2] * u[8] =#
du[1] = 0

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
        n_H * 2e-9 * u[2] * u[13] # added this extra term from a CR ionization reaction
        + 2.0e-10 * Go * exp(-1.9 * Av) * u[14] # this term was added as part of the skipped photoreaction


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
        n_H * 2e-9 * u[2] * u[14] 
        + 2.0e-10 * Go * exp(-1.9 * Av) * u[14] # this term was added as part of the skipped photoreaction

# 14: M
du[14] = n_H * (3.8e-10 * u[13] * u[3]) / (T^0.65) - 
        n_H * 2e-9 * u[2] * u[14] 
        - 2.0e-10 * Go * exp(-1.9 * Av) * u[14] # this term was added as part of the skipped photoreaction

end