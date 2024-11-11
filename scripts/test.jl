#using CUDA, LinearAlgebra
using Revise
using DrWatson
using Catalyst
using DifferentialEquations
using Plots
using OrdinaryDiffEq
using ModelingToolkit
using LSODA


    # %% Set The timespan, paramters, and initial conditions
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs

params = Dict(
    :T => 60,  # from glover paper, section 4
    :n_H => 300, # from glover paper, section 4
    :Av => 1,
    :k1 => 10^(-17.845 + 0.762 * log(T) + 0.1523 * (log(T))^2 - 0.03274 * (log(T))^3),
    :k2 => 0.0000000015,
    :k6 => 0.00000001,
    :k12 => 0.0000000000001269 * ((315614/T)^(1.503)) * (1 + (604625/T)^(0.47))^(-1.923),
    :k14 => 0.0000000025634 * (T)^(1.78186),
    :k15 => 0.0000000069 * T^(-0.35),
    :k17 => 0.0000000001 * T^(-0.5) * (12.72 - 1.615 * log(T) - 0.3162 * (log(T))^2 + 0.0493 * (log(T))^3),
    :k19 => 0.00000000126 * T^(-0.75) * exp(-127500/T),
    :k20 => 0.00000000000467 * (T/300)^(-0.6),
    :k21 => 0.00000000013 * T^(-0.64),
    :k29 => 0.0000000000000000858 * (T)^(0.757),
    :k38 => 0.000000000066,
    :k42 => 0.00000000005 * (T/300)^(0.5),
    :k47 => 0.000000000035,
    :k52 => 0.000000000047 * (T/300)^(-0.34),
    :k146 => 2.1e-19,
    :k149 => 2.5e-18,
    :k154 => 1.32e-32 * (T/300)^-0.38,
    :k157 => 5.99e-33 * (T/5000)^(-1.6), 
    :k158 => 6.16e-29 * (T/300)^(-3.08),
    :R188 => 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R189 => 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R190 => 5.0e-10 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R191 => 1.5e-10 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R192 => 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R193 => 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R194 => 7.5e-12 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R195 => 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2)) 

u0 = [  
    7,        # 1: H⁻ #FIX
    .1,        # 2: H2⁺ #FIX
    9.059e-9, # 3: H3⁺ # This is what nelson has
    .1,        # 4: CH⁺ #FIX
    .1,        # 5: CH2⁺ #FIX
    .1,        # 6: OH⁺ #FIX
    .1,        # 7: H2O⁺ #FIX
    .1,        # 8: H3O⁺ #FIX
    .1,        # 9: CO⁺ #FIX
    0,        # 10: HOC⁺ # This is what nelson has for HCO+
    .1,        # 11: O⁻ #FIX
    .1,        # 12: C⁻ #FIX
    .1,        # 13: O2⁺ #FIX
    2.0e-4,   # 14: e # This is what nelson has
    .1,        # 15: H⁺ #FIX
    .1,        # 16: H #FIX
    0.5,      # 17: H2 # This is what nelson has
    0.1,      # 18: He # This is what nelson has
    7.866e-7, # 19: He⁺ # This is what nelson has
    1.41e-4,  # 20: C # Glover paper has 1.41e-4 section 4, (nelson has 0)
    0.0002,   # 21: C⁺ # This is what nelson has
    3.16e-4,  # 22: O # Glover paper has 3.16e-4 section 4, (nelson has 4e-4)
    .2,        # 23: O⁺ #FIX
    0,        # 24: OH # This is what nelson has, OHx
    0,        # 25: H2O # This is what nelson has, OHxe8
    0,        # 26: CO # This is what nelson has
    .1,        # 27: C2 #FIX
    0,        # 28: O2 # This is what nelson has, OHx
    .1,        # 29: H2O⁺ #FIX
    .1,        # 30: CH #FIX
    0,        # 31: CH2 # This is what nelson has, CHx
    0,        # 32: CH3⁺ # This is what nelson has, CHx
    2.0e-7,   # 33: M # This is what nelson has
    ]


# %% Network admin things
allvars = @strdict tspan u0 params
@variables t 
@species H⁻(t) H2⁺(t) H3⁺(t) CH⁺(t) CH2⁺(t) OH⁺(t) H2O⁺(t) H3O⁺(t) CO⁺(t) HOC⁺(t) O⁻(t) C⁻(t) O2⁺(t) e(t) H⁺(t) H(t) H2(t) He(t) He⁺(t) C(t) C⁺(t) O(t) O⁺(t) OH(t) H2O(t) CO(t) C2(t) O2(t) HCO⁺(t) CH(t) CH2(t) CH3⁺(t) M(t)
#@parameters T n_H

# %% Define the network
reaction_equations = [
    (@reaction n_H * k1, H+e --> H⁻), 
    (@reaction n_H * k2, H⁻+H --> H2+e), 
    (@reaction n_H * 10^(-19.38 - 1.523 * log(T) + 1.118 * (log(T))^2 - 0.1269 * (log(T))^3) , H+ H⁺--> H2⁺ ), 
    (@reaction n_H * 0.00000000064 , H+H2⁺ --> H2 + H⁺), 
    (@reaction n_H * 0.0000024 * (T)^(-1/2) * (1 + T/20000) , H⁻ + H⁺--> H+H), 
    (@reaction n_H * k6, H2+ + e --> H+H), 
    (@reaction n_H * (-0.00000033232183 + 0.00000033735382 * log(T) - 0.00000014491368 * (log(T))^2 + 0.000000034172805 * (log(T))^3 - 0.000000004781372 * (log(T))^4 + 0.00000000039731542 * (log(T))^5 - 0.000000000018171411 * (log(T))^6 + 35311932000000 * (log(T))^7) * exp(-21237.15/T) , H2+H⁺--> H2⁺+H), 
    (@reaction n_H * 0.00000000373 * (T)^(0.1121) * exp(-99430/T) , H2+e --> H+H+e), 
    (@reaction n_H * 0.00000000000667 * (T)^(1/2) * exp(-(1 + 63590/T)) , H2+H --> H+H+H), 
    (@reaction n_H * (5.996E+30 * (T)^(4.1881)) / ((1 + 0.000006761 * T)^(5.6881)) * exp(-54657.4/T) , H2+H2 --> H2+H+H), 
    (@reaction n_H * exp(-32.71396786 + 13.536556 * log(T) - 5.73932875 * (log(T))^2 + 1.56315498 * (log(T))^3 - 0.2877056 * (log(T))^4 + 0.0348255977 * (log(T))^5 - 0.00263197617 * (log(T))^6 + 0.000111954395 * (log(T))^7 - 0.00000203914985 * (log(T))^8) , H+e --> H⁺ + e+e), 
    (@reaction n_H * k12, H⁺ + e --> H ),
    (@reaction n_H * exp(-18.01849334 + 2.3608522 * log(T) - 0.2827443 * (log(T))^2 + 0.0162331664 * (log(T))^3 - 0.0336501203 * (log(T))^4 + 0.0117832978 * (log(T))^5 - 0.0016561947 * (log(T))^6 + 0.00010682752 * (log(T))^7 - 0.00000263128581 * (log(T))^8) , H⁻ + e --> H + e + e), 
    (@reaction n_H * k14, H⁻ + H --> H + H + e), 
    (@reaction n_H * k15, H⁻ + H⁺--> H2⁺+e),
    (@reaction n_H * exp(-44.09864886 + 23.91596563 * log(T) - 10.7532302 * (log(T))^2 + 3.05803875 * (log(T))^3 - 0.56851189 * (log(T))^4 + 0.0679539123 * (log(T))^5 - 0.0050090561 * (log(T))^6 + 0.000206723616 * (log(T))^7 - 0.00000364916141 * (log(T))^8) , He+e --> He⁺ + e+e), 
    (@reaction n_H * k17, He⁺ + e --> He ), 
    (@reaction n_H * 0.00000000000000125 * (T/300)^(0.25) , He⁺ + H --> He+H⁺), 
    (@reaction n_H * k19, He+H⁺--> He⁺ + H), 
    (@reaction n_H * k20, C⁺ + e --> C ), 
    (@reaction n_H * k21, O⁺ + e --> O ), 
    (@reaction n_H * 0.0000000685 * (0.193 + (11.26/T))^(-1) * (11.26/T)^(0.25) * exp(-(11.26/T)) , C+e --> C⁺ + e+e), 
    (@reaction n_H * 0.0000000359 * (0.073 + (13.6/T))^(-1) * (13.6/T)^(0.34) * exp(-13.6/T) , O+e --> O⁺ + e+e), 
    (@reaction n_H * 0.0000000000499 * T^(0.405) + 0.000000000754 * T^(-0.458) , O⁺ + H --> O+H⁺), 
    (@reaction n_H * (0.0000000000108 * T^(0.517) + 0.0000000004 * T^(0.00669)) * exp(-227/T) , O+H⁺--> O⁺ + H), 
    (@reaction n_H * 0.000000000000004991 * (T/10000)^(0.3794) * exp(-T/1121000) + 0.00000000000000278 * (T/10000)^(-0.2163) * exp(T/815800) , O+He⁺ --> O⁺ + He), 
    (@reaction n_H * 0.00000000000000039 * T^(0.213) , C+H⁺--> C⁺ + H), 
    (@reaction n_H * 0.0000000000000608 * (T/10000)^(1.96) * exp(-170000/T) , C⁺ + H --> C+H⁺), 
    (@reaction n_H * k29, C+He⁺ --> C⁺ + He), 
    (@reaction n_H * 10^(-27.029 + 3.801 * log(T) - 29487/T) , H2+He --> H+H+He), 
    (@reaction n_H * 0.000000006 * exp(-50900/T) , OH+H --> O+H+H), 
    (@reaction n_H * 0.00000000038 , HOC⁺ + H2 --> HCO⁺ + H2), 
    (@reaction n_H * 0.0000000004 , HOC⁺ + CO --> HCO⁺ + CO), 
    (@reaction n_H * 0.000000000664 * exp(-11700/T) , C+H2 --> CH+H), 
    (@reaction n_H * 0.000000000131 * exp(-80/T) , CH+H --> C+H2), 
    (@reaction n_H * 0.000000000546 * exp(-1943/T) , CH+H2 --> CH2+H), 
    (@reaction n_H * 0.0000000000659 , CH+C --> C2+H), 
    (@reaction n_H * k38, CH+O --> CO+H), 
    (@reaction n_H * 0.0000000000664 , CH2+H --> CH+H2), 
    (@reaction n_H * 0.000000000133 , CH2+O --> CO+H+H), 
    (@reaction n_H * 0.00000000008 , CH2+O --> CO+H2), 
    (@reaction n_H * k42, C2+O --> CO+C), 
    (@reaction n_H * 0.000000000000314 * (T/300)^(2.7) * exp(-3150/T) , O+H2 --> OH+H), 
    (@reaction n_H * 0.0000000000000699 * (T/300)^(2.8) * exp(-1950/T) , OH+H --> O+H2), 
    (@reaction n_H * 0.00000000000205 * (T/300)^(1.52) * exp(-1736/T) , OH+H2 --> H2O+H), 
    (@reaction n_H * 0.0000000001 , OH+C --> CO+H), 
    (@reaction n_H * k47, OH+O --> O2+H), 
    (@reaction n_H * 0.00000000000165 * (T/300)^(1.14) * exp(-50/T) , OH+OH --> H2O+H), 
    (@reaction n_H * 0.0000000000159 * (T/300)^(1.2) * exp(-9610/T) , H2O+H --> H2+OH), 
    (@reaction n_H * 0.000000000261 * exp(-8156/T) , O2+H --> OH+O), 
    (@reaction n_H * 0.000000000316 * exp(-21890/T) , O2+H2 --> OH+OH), 
    (@reaction n_H * k52, O2+C --> CO+O), 
    (@reaction n_H * 0.00000000011 * (T/300) * exp(-77700/T) , CO+H --> C+OH), 
    (@reaction n_H * 0.00000000224 * (T/300)^(0.042) * exp(-T/46600) , H2⁺+H2 --> H3⁺+H), 
    (@reaction n_H * 0.0000000077 * exp(-17560/T) , H3⁺+H --> H2⁺+H2), 
    (@reaction n_H * 0.0000000024 , C+H2⁺ --> CH⁺ + H), 
    (@reaction n_H * 0.000000002 , C+H3⁺ --> CH⁺ + H2), 
    (@reaction n_H * 0.0000000001 * exp(-4640/T) , C⁺ + H2 --> CH⁺ + H),
    (@reaction n_H * 0.00000000075 , CH⁺ + H --> C⁺ + H2), 
    (@reaction n_H * 0.0000000012 , CH⁺ + H2 --> CH2⁺+H), 
    (@reaction n_H * 0.00000000035 , CH⁺ + O --> CO⁺ + H), 
    (@reaction n_H * 0.0000000014 , CH2+H⁺--> CH⁺ + H2), 
    (@reaction n_H * 0.000000001 * exp(-7080/T) , CH2⁺+H --> CH⁺ + H2), 
    (@reaction n_H * 0.0000000016 , CH2⁺+H2 --> CH3⁺+H), 
    (@reaction n_H * 0.00000000075 , CH2⁺+O --> HCO⁺ + H), 
    (@reaction n_H * 0.0000000007 * exp(-10560/T) , CH3⁺+H --> CH2⁺+H2), 
    (@reaction n_H * 0.0000000004 , CH3⁺+O --> HCO⁺ + H2), 
    (@reaction n_H * 0.00000000048 , C2+O⁺ --> CO⁺ + C), 
    (@reaction n_H * 0.0000000017 , O⁺ + H2 --> OH⁺ + H), 
    (@reaction n_H * 0.0000000015 , O+H2⁺ --> OH⁺ + H), 
    (@reaction n_H * 0.00000000084 , O+H3⁺ --> OH⁺ + H2), 
    (@reaction n_H * 0.0000000013 , OH+H3⁺ --> H2O⁺ + H2), 
    (@reaction n_H * 0.00000000077 , OH+C⁺ --> CO⁺ + H), 
    (@reaction n_H * 0.00000000101 , OH⁺ + H2 --> H2O⁺ + H), 
    (@reaction n_H * 0.00000000064 , H2O⁺ + H2 --> H3O⁺ + H), 
    (@reaction n_H * 0.0000000059 , H2O+H3⁺ --> H3O⁺ + H2), 
    (@reaction n_H * 0.0000000009 , H2O+C⁺ --> HCO⁺ + H), 
    (@reaction n_H * 0.0000000018 , H2O+C⁺ --> HOC⁺ + H), 
    (@reaction n_H * 0.00000000001 , H3O⁺ + C --> HCO⁺ + H2), 
    (@reaction n_H * 0.00000000038 , O2+C⁺ --> CO⁺ + O), 
    (@reaction n_H * 0.00000000062 , O2+C⁺ --> CO+O⁺), 
    (@reaction n_H * 0.00000000091 , O2+CH2⁺ --> HCO⁺ + OH), 
    (@reaction n_H * 0.000000000052 , O2⁺+C --> CO⁺ + O), 
    (@reaction n_H * 0.000000000027 , CO+H3⁺ --> HOC⁺ + H2), 
    (@reaction n_H * 0.0000000017 , CO+H3⁺ --> HCO⁺ + H2), 
    (@reaction n_H * 0.0000000011 , HCO⁺ + C --> CO+CH⁺), 
    (@reaction n_H * 0.0000000025 , HCO⁺ + H2O --> CO+H3O⁺), 
    (@reaction n_H * 0.0000000000000072 , H2+He⁺  -->  He+H2⁺), 
    (@reaction n_H * 0.000000000000037 * exp(-35/T) , H2+He⁺  -->  He+H+H⁺), 
    (@reaction n_H * 0.0000000019 , CH+H⁺ -->  CH⁺ + H), 
    (@reaction n_H * 0.0000000014 , CH2+H⁺ -->  CH2⁺+H), 
    (@reaction n_H * 0.00000000075 , CH2+He⁺ -->  C⁺ + He+H2), 
    (@reaction n_H * 0.0000000016 , C2+He⁺ -->  C⁺ + C+He), 
    (@reaction n_H * 0.0000000021 , OH+H⁺ -->  OH⁺ + H), 
    (@reaction n_H * 0.0000000011 , OH+He⁺ -->  O⁺ + He+H), 
    (@reaction n_H * 0.0000000069 , H2O+H⁺ -->  H2O⁺ + H), 
    (@reaction n_H * 0.000000000204 , H2O+He⁺ -->  OH+He+H⁺), 
    (@reaction n_H * 0.000000000286 , H2O+He⁺ -->  OH⁺ + He+H), 
    (@reaction n_H * 0.0000000000605 , H2O+He⁺ -->  H2O⁺ + He), 
    (@reaction n_H * 0.000000002 , O2+H⁺ -->  O2⁺+H), 
    (@reaction n_H * 0.000000000033 , O2+He⁺ -->  O2⁺+He), 
    (@reaction n_H * 0.0000000011 , O2+He⁺ -->  O⁺ + O+He), 
    (@reaction n_H * 0.000000000052 , O2⁺+C  -->  O2+C⁺), 
    (@reaction n_H * 10000-9*(T/300)^(-0.5) , CO+He⁺ -->  C⁺ + O+He), 
    (@reaction n_H * 0.00000000000000014 * (T/300)^(-0.5) , CO+He⁺ -->  C+O⁺ + He), 
    (@reaction n_H * 0.00000000075 , CO⁺ + H  -->  CO+H⁺), 
    (@reaction n_H * 0.00000023 * (T/300)^(-0.5) , C⁻ + H⁺ -->  C + H), 
    (@reaction n_H * 0.00000023 * (T/300)^(-0.5) , O⁻ + H⁺ -->  O + H), 
    (@reaction n_H * 0.000000232 * (T/300)^(-0.52) * exp(T/22400) , He⁺ + H⁻  -->  He+H), 
    (@reaction n_H * 0.0000000234 * (T/300)^(-0.52) , H3⁺+e  -->  H2+H), 
    (@reaction n_H * 0.0000000436 * (T/300)^(-0.52) , H3⁺+e  -->  H+H+H), 
    (@reaction n_H * 0.00000007 * (T/300)^(-0.5) , CH⁺ + e  -->  C+H), 
    (@reaction n_H * 0.00000016 * (T/300)^(-0.6) , CH2⁺+e  -->  CH+H), 
    (@reaction n_H * 0.000000403 * (T/300)^(-0.6) , CH2⁺+e  -->  C+H+H), 
    (@reaction n_H * 0.0000000768 * (T/300)^(-0.6) , CH2⁺+e  -->  C+H2), 
    (@reaction n_H * 0.0000000775 * (T/300)^(-0.5) , CH3⁺+e  -->  CH2+H), 
    (@reaction n_H * 0.000000195 * (T/300)^(-0.5) , CH3⁺+e  -->  CH+H2), 
    (@reaction n_H * 0.0000002 * (T/300)^(-0.4) , CH3⁺+e  -->  CH+H+H), 
    (@reaction n_H * 0.0000000063 * (T/300)^(-0.48) , OH⁺ + e  -->  O+H), 
    (@reaction n_H * 0.000000305 * (T/300)^(-0.5) , H2O⁺ + e  -->  O+H+H), 
    (@reaction n_H * 0.000000039 * (T/300)^(-0.5) , H2O⁺ + e  -->  O+H2), 
    (@reaction n_H * 0.000000086 * (T/300)^(-0.5) , H2O⁺ + e  -->  OH+H), 
    (@reaction n_H * 0.000000108 * (T/300)^(-0.5) , H3O⁺ + e  -->  H+H2O), 
    (@reaction n_H * 0.0000000602 * (T/300)^(-0.5) , H3O⁺ + e  -->  OH+H2), 
    (@reaction n_H * 0.000000258 * (T/300)^(-0.5) , H3O⁺ + e  -->  OH+H+H), 
    (@reaction n_H * 0.0000000056 * (T/300)^(-0.5) , H3O⁺ + e  -->  O+H+H2), 
    (@reaction n_H * 0.000000195 * (T/300)^(-0.7) , O2⁺+e  -->  O+O), 
    (@reaction n_H * 0.000000275 * (T/300)^(-0.55) , CO⁺ + e  -->  C+O), 
    (@reaction n_H * 0.000000276 * (T/300)^(-0.64) , HCO⁺ + e  -->  CO+H), 
    (@reaction n_H * 0.000000024 * (T/300)^(-0.64) , HCO⁺ + e  -->  OH+C), 
    (@reaction n_H * 0.00000011 * (T/300)^(-1) , HOC⁺ + e  -->  CO+H), 
    (@reaction n_H * 0.000000001 , H⁻+C  -->  CH+e), 
    (@reaction n_H * 0.000000001 , H⁻+O  -->  OH+e), 
    (@reaction n_H * 0.0000000001 , H⁻+OH  -->  H2O+e), 
    (@reaction n_H * 0.0000000005 , C⁻ + H  -->  CH + e), 
    (@reaction n_H * 0.0000000000001 , C⁻ + H2  -->  CH2 + e), 
    (@reaction n_H * 0.0000000005 , C⁻ + O  -->  CO + e), 
    (@reaction n_H * 0.0000000005 , O⁻ + H  -->  OH + e), 
    (@reaction n_H * 0.0000000007 , O⁻ + H2  -->  H2O + e), 
    (@reaction n_H * 0.0000000005 , O⁻ + C  -->  CO + e), 
    (@reaction n_H * 0.0000000000000001 , H2+ H⁺ -->  H3⁺ ), 
    (@reaction n_H * 0.00000000000000225 , C+e  -->  C⁻ ), 
    (@reaction n_H * 0.00000000000000001 , C+H  -->  CH ), 
    (@reaction n_H * 0.00000000000000001 , C+H2  -->  CH2 ), 
    (@reaction n_H * 4360000000000000000 * (T/300)^(0.35) * exp(-161.3/T) , C+C  -->  C2 ), 
    (@reaction n_H * k146, C + O  --> CO ), 
    (@reaction n_H * 0.000000000000000446 * T^(-0.5) * exp(-4.93 / (T^(2/3)) ) , C⁺ + H  -->  CH⁺  ), 
    (@reaction n_H * 0.0000000000000004 *(T/300)^(-0.2) , C⁺ + H2  -->  CH2⁺ ), 
    (@reaction n_H * k149, C⁺ + O  -->  CO⁺  ), 
    (@reaction n_H * 0.0000000000000015 , O+e  -->  O⁻ ), 
    (@reaction n_H * 9.9e-19*(T/300)^(-0.38) , O+H  -->  OH ), 
    (@reaction n_H * 4.9e-20 * (T/300)^1.58 , O+O  -->  O2 ), 
    (@reaction n_H * 5.26e-18 *(T/300)^(-5.22) *exp(-90/T) , OH+H  -->  H2O ), 
    (@reaction n_H * k154, H+H+H  -->  H2+H), 
    (@reaction n_H * 2.8e-31 * (T)^-0.6 , H+H+H2  -->  H2+H2), 
    (@reaction n_H * 6.9e-32 * (T)^(-0.4) , H+H+He  -->  H2+He), 
    (@reaction n_H * k157, C+C+M  -->  C2+M), 
    (@reaction n_H * k158, C+O+M  -->  CO+M), 
    (@reaction n_H * 4.33e-32 * (T/300)^(-1) , O+H+M  -->  OH+M), 
    (@reaction n_H * 2.56e-31 * (T/300)^(-2) , OH+H+M  -->  H2O+M), 
    (@reaction n_H * 9.2e-34 * (T/300)^(-1) , O+O+M  -->  O2+M), 
    (@reaction n_H * 0.00000000002 * (T/300)^(0.44) , O+CH  -->  HCO⁺ + e), 
    (@reaction 7.1e-7 * exp(-0.5 * Av), H⁻    -->  H+e), 
    (@reaction 1.1e-9 * exp(-1.9 * Av), H2⁺   -->  H+H⁺), 
    (@reaction 4.9e-13 * exp(-1.8 * Av), H3⁺   -->  H2+H⁺),
    (@reaction 4.9e-13 * exp(-2.3 * Av), H3⁺   -->  H2⁺+H), 
    (@reaction 3.1e-10 * exp(-3.0 * Av), C   -->  C⁺ + e), 
    (@reaction 2.4e-7 * exp(-0.9 * Av), C⁻  -->  C+e), 
    (@reaction 8.7e-10 * exp(-1.2 * Av), CH   -->  C+H), 
    (@reaction 7.7e-10 * exp(-2.8 * Av), CH   -->  CH⁺ + e), 
    (@reaction 2.6e-10 * exp(-2.5 * Av), CH⁺    -->  C+H⁺), 
    (@reaction 7.1e-10 * exp(-1.7 * Av), CH2   -->  CH+H), 
    (@reaction 5.9e-10 * exp(-2.3 * Av), CH2   -->  CH2⁺+e), 
    (@reaction 4.6e-10 * exp(-1.7 * Av), CH2⁺   -->  CH⁺ + H), 
    (@reaction 1.0e-9 * exp(-1.7 * Av), CH3⁺   -->  CH2⁺+H), 
    (@reaction 1.0e-9 * exp(-1.7 * Av), CH3⁺   -->  CH⁺ + H2), 
    (@reaction 1.5e-10 * exp(-2.1 * Av), C2   -->  C+C), 
    (@reaction 2.4e-7 * exp(-0.5 * Av), O⁻ -->  O+e),
    (@reaction 3.7e-10 * exp(-1.7 * Av), OH   -->  O+H), 
    (@reaction 1.6e-12 * exp(-3.1 * Av), OH   -->  OH⁺ + e), 
    (@reaction 1.0e-12 * exp(-1.8 * Av), OH⁺    -->  O+H⁺), 
    (@reaction 6.0e-10 * exp(-1.7 * Av), H2O   -->  OH+H), 
    (@reaction 3.2e-11 * exp(-3.9 * Av), H2O   -->  H2O⁺ + e), 
    (@reaction R188, H2O⁺    -->  H2⁺+O), 
    (@reaction R189, H2O⁺    -->  H⁺ + OH), 
    (@reaction R190, H2O⁺    -->  O⁺ + H2), 
    (@reaction R191, H2O⁺    -->  OH⁺ + H), 
    (@reaction R192, H3O⁺    -->  H⁺ + H2O), 
    (@reaction R193, H3O⁺    -->  H2⁺+OH), 
    (@reaction R194, H3O⁺    -->  H2O⁺ + H), 
    (@reaction R195, H3O⁺    -->  OH⁺ + H2), 
    (@reaction 5.6e-11 * exp(-3.7 * Av), O2   -->  O2⁺+e), 
    (@reaction 7.0e-10 * exp(-1.8 * Av), O2   -->  O+O), 
    (@reaction 1 , H   -->  H⁺ + e), 
    (@reaction 1.1 , He   -->  He⁺ + e), 
    (@reaction 0.037 , H2   -->  H⁺ + H+e), 
    (@reaction 0.22 , H2   -->  H+H), 
    (@reaction 0.00065 , H2   -->  H⁺), 
    (@reaction 2 , H2   -->  H2⁺+e), 
    (@reaction 3.8 , C   -->  C⁺ + e), 
    (@reaction 5.7 , O   -->  O⁺ + e), 
    (@reaction 6.5 , CO   -->  CO⁺ + e), 
    (@reaction 2800 , C   -->  C⁺ + e), 
    (@reaction 4000 , CH   -->  C+H), 
    (@reaction 960 , CH⁺    -->  C⁺ + H), 
    (@reaction 2700 , CH2   -->  CH2⁺+e), 
    (@reaction 2700 , CH2   -->  CH+H), 
    (@reaction 1300 , C2   -->  C+C), 
    (@reaction 2800 , OH   -->  O+H), 
    (@reaction 5300 , H2O   -->  OH+H), 
    (@reaction 4100 , O2   -->  O+O), 
    (@reaction 640  , O2   -->  O2⁺ + e), 
    
]


    # %% Turn the network into an ODE system
@named system = ReactionSystem(reaction_equations, t)
odesys = convert(ODESystem, complete(system))
sys = convert(ODESystem, complete(system))
ssys = structural_simplify(sys)

prob = ODEProblem(ssys, u0, tspan, params)
#sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
sol = solve(prob, Rodas4())
plot(sol, idxs = (0,12), lw = 3, title = "HCO+ from Glover network")
#14 is HOC+
