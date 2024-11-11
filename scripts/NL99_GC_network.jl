    # %% Set The timespan, paramters, and initial conditions
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs

params = Dict(
    :T => 10,
    :n_H => 611,
    :Go => 1.7,
    :Av => 2, 
    :shield_CO => 1, 
    :shield_H2 => 1)

u0 = [0.5,      # 1:  H2   yep?
        1,        # 2:  H    fix
        1,        # 3:  H+   fix
        9.059e-9, # 4:  H3+  yep
        2.0e-4,   # 5:  e    yep
        0.1,      # 6:  He  SEE lines 535 NL99
        7.866e-7, # 7:  He+  yep? should be 2.622e-5
        0.0,      # 8:  C    yep
        0.0,      # 9:  CHx  yep
        0.0004,   # 10:  O    yep
        0.0,      # 11:  OHx  yep
        0.0,      # 12: CO   yep
        0.0,      # 13: HCO+ yep
        0.0002,   # 14: C+   yep
        2.0e-7,   # 15: M+   yep
        2.0e-7]   # 16: M    yep





    # %% Network admin things
allvars = @strdict params tspan u0
@variables t 
@species H2(t) H(t) H⁺(t) H3⁺(t) el(t) He(t) He⁺(t) C(t) CHx(t) O(t) OHx(t) CO(t) HCO⁺(t) C⁺(t) M⁺(t) M(t)
@parameters T n_H Go Av shield_CO shield_H2

    # %% Define the network
reaction_equations = [
    # Cosmic-ray Ionization
    (@reaction 6.0e-18, H --> H⁺ + el),
    (@reaction 1.2e-17, H2 --> H3⁺ + H + el),
    (@reaction 6.8e-18, He --> He⁺ + el),

    # Photoreactions
    (@reaction 3.3e-11 * shield_H2 * Go * exp(-3.74 * Av), H2 --> H + H),  # added reaction
    (@reaction 1.2e-10 * shield_CO * Go * exp(-3.53 * Av), CO --> C + O), #updated rate: 1e-10 -> 1.2e-10
    (@reaction 1.8e-10 * Go * exp(-3 * Av), C --> C⁺ + el), # updated rate: 3e-10 -> 1.8e-10
    (@reaction 1e-9 * Go * exp(-1.5 * Av), CHx --> C + H), 
    (@reaction 5e-10 * Go * exp(-1.7 * Av), OHx --> O + H), 
    (@reaction 2e-10 * Go * exp(-1.9 * Av), M --> M⁺ + el), 
    (@reaction 1.5e-10 * Go * exp(-2.5 * Av), HCO⁺ --> CO + H⁺),  # is nelson wrong? Nelson days its HCO+ --> CO + H

    # Two Body Reactions
    (@reaction n_H * 2e-9, H3⁺ + C --> CHx + H2), 
    (@reaction n_H * 8e-10, H3⁺  + O   --> OHx + H2),
    (@reaction n_H * 1.7e-9, H3⁺  + CO  --> HCO⁺ + H2),
    (@reaction n_H * 7e-15, He⁺  + H2  --> He + H + H⁺),
    (@reaction n_H * 1.4e-9*((T/300)^(-0.5)), He⁺  + CO  --> C⁺ + O + He),  #Update rate: 1.6e-9 -> 1.4e-9*(T/300.)**-0.5
    (@reaction n_H * 4e-16, C⁺ + H2  --> CHx + H),
    (@reaction n_H * 1e-9, C⁺ + OHx --> HCO⁺),
    (@reaction n_H * 2e-10, O + CHx --> CO + H),
    (@reaction n_H * 5.8e-12 * (T)^(0.5), C + OHx --> CO + H),
    (@reaction n_H * (1.0e-11/sqrt(T)) * (11.19 - 1.676*log(T) - 0.2852*log(T)^2 + 0.04433*log(T)^3), He⁺ + el  --> He), #Updated rate: 9e-11 * (T)^(-0.64) -> 1.0e-11/np.sqrt(T) * \ (11.19 - 1.676*logT - 0.2852*logT**2 + 0.04433*logT**3)
    (@reaction n_H * 2.34e-8 * (T/300)^(-0.52), H3⁺ + el  --> H2 + H), # Updated rate: 1.9e-6 * (T)^(-.54) ->2.34e-8 * (T/300.)**-0.52
    (@reaction n_H * 4.36e-8 * (T/300)^(-0.52), H3⁺ + el  --> 3H), # new reaction
    (@reaction n_H * 4.67e-12 * (T/300)^(-0.6), C⁺ + el  --> C), # Updated rate: 1.4e-10 * (T)^(-0.61) -> 4.67e-12 * (T/300.)**-0.6
    (@reaction n_H * 2.76e-7 * (T/300)^(-0.64), HCO⁺ + el  --> CO + H), # updated rate: 3.3e-5 * (T)^(-1.0) -> 2.76e-7 * (T/300.)**-0.64
    (@reaction n_H * 2.753e-14*(315614/T)^(1.5) * (1.0+(115188/T)^(0.407))^(-2.242), H⁺ + el  --> H), # New reaction
    #15 below FIX
    (@reaction n_H * 1, H2 + H   --> 3H), # New reaction
    (@reaction n_H * 1, H2 + H2  --> H2 + 2H), # New reaction
    (@reaction n_H * 1, H + el  --> H⁺ + 2el), # New reaction
    (@reaction n_H * 1, He⁺ + H2  --> H⁺ + He + H), # New reaction
    (@reaction n_H * 1, H + H   --> H2), # Be careful with grain = true
    (@reaction n_H * 1, H⁺ + el  --> H), # Be careful with grain = true    
    # finish here
    (@reaction n_H * 3.8e-10 * (T)^(-0.65), M⁺ + el  --> M),
    (@reaction n_H * 2e-9, H3⁺ + M   --> M⁺ + H + H2),
    ]


    # %% Turn the network into an ODE system
@named system = ReactionSystem(reaction_equations, t)
odesys = convert(ODESystem, complete(system))
sys = convert(ODESystem, complete(system))
ssys = structural_simplify(sys)