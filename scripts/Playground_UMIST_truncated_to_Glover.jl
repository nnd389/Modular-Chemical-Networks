### Packages ###
using Catalyst
using DifferentialEquations
using Plots
using Sundials, LSODA
using OrdinaryDiffEq
using Symbolics
using DiffEqDevTools
using ODEInterface, ODEInterfaceDiffEq
using ModelingToolkit


### Timespan ###
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs

### Initial Conditions ###
u0 = [  0,        # 1: H
        2e-4,     # 2: e Nelson has 2e-4
        0,        # 3: H-
        0.5,      # 4: H2 Nelson has 0.5
        0,        # 5: H+
        0,        # 6: H2+
        0.1,      # 7: He Nelson has 0.1
        7.866e-7, # 8: He+ Nelson has 7.866e-7
        2e-4,     # 9: C+ Nelson has 2e-4, Glover has...
        0,        # 10: C Nelson has 0, Glover has 1.41e-4, see section 4
        0,        # 11: O+
        4e-4,     # 12: O Nelson has 4e-4, Glover has 3.16e-4, see section 4
        0,        # 13: OH Nelson has 0 for OHx
        0,        # 14: HOC+
        0,        # 15: HCO+ Nelson has 0
        0,        # 16: CO Nelson has 0
        0,        # 17: CH
        0,        # 18: CH2 Nelson has 0 for CHx
        0,        # 19: C2
        0,        # 20: H2O Nelson has 0
        0,        # 21: O2 Nelson has 0, OHx FIX-------CHECK IF THIS IS O2+ FIND O2
        9.059e-9, # 22: H3+ Nelson has 9.059e-9
        0,        # 23: CH+
        0,        # 24: CH2+
        0,        # 25: CO+
        0,        # 26: CH3+ Nelson has 0 for CHx
        0,        # 27: OH+
        0,        # 28: H2O+
        0,        # 29: H3O+
        0,        # 30: O2+
        0,        # 31: C-
        0,        # 32: O-
        2e-7]     # 33: M Nelson has 2e-7
        

### Parameters ### Check 15 and 48 and 105 and 130

### Parameters ###
# UMIST is using 5.98e-18 as their cosmic ionization rate
# Glover uses omega = 0.6
T = 10
omega = 0.5
Av = 2
n_H = 611
cr_ion_rate = 6e-18

params = Dict(
         :T => T,
         :Av => Av,
         :n_H => n_H,
         :cr_ion_rate => cr_ion_rate,
         :omega => omega,
         :k1 =>  n_H * 3.37e-16 * (T/300)^(0.64) * exp(-(9.2/T)),     # Found Match for Glover 1 H + e --> H⁻ and Umist 6155 H + e --> H⁻ with temp range [10, 41000]
         :k2 =>  n_H * 4.82e-09 * (T/300)^(0.02) * exp(-(4.3/T)),     # Found Match for Glover 2 H⁻ + H --> H2 + e and Umist 75 H⁻ + H --> H2 + e with temp range [10, 100]
         :k3 =>  n_H * 1.15e-18 * (T/300)^(1.49) * exp(-(228.0/T)),     # Found Match for Glover 3 H + H⁺ --> H2⁺  and Umist 6093 H⁺ + H --> H2⁺ with temp range [200, 32000]
         :k4 =>  n_H * 6.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 4 H + H2⁺ --> H2 + H⁺ and Umist 489 H + H2⁺ --> H2 + H⁺ with temp range [10, 41000]
         :k5 =>  n_H * 7.51e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 5 H⁻ + H⁺--> H + H and Umist 4921 H⁻ + H⁺ --> H + H with temp range [10, 41000]
         :k6 =>  n_H * 1.6e-08 * (T/300)^(-0.43) * exp(-(0.0/T)),     # Found Match for Glover 6 H2⁺ + e --> H + H and Umist 1239 H2⁺ + e --> H + H with temp range [10, 300]
         # :k7 =>  Glover 7 H2 + H⁺ --> H2⁺ + H had 0 product matches
         :k8 =>  n_H * 3.22e-09 * (T/300)^(0.35) * exp(-(102000.0/T)),     # Found Match for Glover 8 H2 + e --> H + H + e and Umist 140 H2 + e --> H + H + e with temp range [3400, 41000]
         :k9 =>  n_H * 4.67e-07 * (T/300)^(-1.0) * exp(-(55000.0/T)),     # Found Match for Glover 9 H2 + H --> H + H + H and Umist 142 H + H2 --> H + H + H with temp range [1833, 41000]
         :k10 =>  n_H * 1e-08 * (T/300)^(0.0) * exp(-(84100.0/T)),     # Found Match for Glover 10 H2 + H2 --> H2 + H + H and Umist 135 H2 + H2 --> H2 + H + H with temp range [2803, 41000]
         # :k11 => NO match for Glover reaction 11 H + e --> H⁺ + e + e
         :k12 =>  n_H * 3.5e-12 * (T/300)^(-0.75) * exp(-(0.0/T)),     # Found Match for Glover 12 H⁺ + e --> H  and Umist 6162 H⁺ + e --> H with temp range [10, 20000]
         # :k13 =>  Glover 13 H⁻ + e --> H + e + e had 0 product matches
         # :k14 => NO match for Glover reaction 14 H⁻ + H --> H + H + e
         # :k15 => FIX
         # :k16 =>  Glover 16 He + e --> He⁺ + e + e had 0 product matches
         :k17 =>  n_H * 5.36e-12 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 17 He⁺ + e --> He  and Umist 6166 He⁺ + e --> He with temp range [10, 1000]
         :k18 =>  n_H * 1.2e-15 * (T/300)^(0.25) * exp(-(0.0/T)),     # Found Match for Glover 18 He⁺ + H --> He + H⁺ and Umist 491 H + He⁺ --> He + H⁺ with temp range [10, 41000]
         # :k19 => NO match for Glover reaction 19 He + H⁺--> He⁺ + H
         


         #k20 => n_H * 4.67e-12 * (T/300)^(-0.6), C⁺ + e --> C ), # Glover rate for C⁺ + e --> C
         :k20 =>  n_H * 2.36e-12 * (T/300)^(-0.29) * exp(-(-17.6/T)),     # Umist rate for C⁺ + e --> C with temp range [10, 41000]
         #:k20 => n_H * 1.4e-10 * (T)^(-0.61), # Nelson rate 1 for C⁺ + e --> C

         :k20_backwards => 2.3e-17, # Umist rate for C --> C⁺ + e 
         :k20_backwards => 3e-10 * 2 * exp(-3 * Av), # Nelson rate 2 for C --> C⁺ + e 





         
         :k21 =>  n_H * 3.24e-12 * (T/300)^(-0.66) * exp(-(0.0/T)),     # Found Match for Glover 21 O⁺ + e --> O  and Umist 6170 O⁺ + e --> O with temp range [10, 41000]
         # :k22 => NO match for Glover reaction 22 C + e --> C⁺ + e + e
         # :k23 => NO match for Glover reaction 23 O + e --> O⁺ + e + e
         :k24 =>  n_H * 5.66e-10 * (T/300)^(0.36) * exp(-(-8.6/T)),     # Found Match for Glover 24 O⁺ + H --> O + H⁺ and Umist 492 H + O⁺ --> O + H⁺ with temp range [10, 41000]
         :k25 =>  n_H * 6.86e-10 * (T/300)^(0.26) * exp(-(224.3/T)),     # Found Match for Glover 25 O + H⁺--> O⁺ + H and Umist 406 H⁺ + O --> O⁺ + H with temp range [10, 41000]
         # :k26 =>  Glover 26 O + He⁺ --> O⁺ + He had 0 product matches
         # :k27 =>  Glover 27 C + H⁺ --> C⁺ + H had 0 product matches
         # :k28 => NO match for Glover reaction 28 C⁺ + H --> C + H⁺
         :k29 =>  n_H * 6.3e-15 * (T/300)^(0.75) * exp(-(0.0/T)),     # Found Match for Glover 29 C + He⁺ --> C⁺ + He and Umist 513 He⁺ + C --> C⁺ + He with temp range [10, 300]
         # :k30 =>  Glover 30 H2 + He --> H + H + He had 0 product matches
         :k31 =>  n_H * 6e-09 * (T/300)^(0.0) * exp(-(50900.0/T)),     # Found Match for Glover 31 OH + H --> O + H + H and Umist 145 H + OH --> O + H + H with temp range [1696, 41000]
         :k32 =>  n_H * 3.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 32 HOC⁺ + H2 --> HCO⁺ + H2 and Umist 137 H2 + HOC⁺ --> HCO⁺ + H2 with temp range [20, 41000]
         # :k33 =>  Glover 33 HOC⁺ + CO --> HCO⁺ + CO had 0 product matches
         :k34 =>  n_H * 6.64e-10 * (T/300)^(0.0) * exp(-(11700.0/T)),     # Found Match for Glover 34 C + H2 --> CH + H and Umist 5348 H2 + C --> CH + H with temp range [300, 2500]
         :k35 =>  n_H * 1.31e-10 * (T/300)^(0.0) * exp(-(80.0/T)),     # Found Match for Glover 35 CH + H --> C + H2 and Umist 5374 H + CH --> C + H2 with temp range [300, 2000]
         :k36 =>  n_H * 5.46e-10 * (T/300)^(0.0) * exp(-(1943.0/T)),     # Found Match for Glover 36 CH + H2 --> CH2 + H and Umist 5351 H2 + CH --> CH2 + H with temp range [300, 2500]
         :k37 =>  n_H * 6.59e-11 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 37 CH + C --> C2 + H and Umist 5173 C + CH --> C2 + H with temp range [10, 300]
         :k38 =>  n_H * 6.02e-11 * (T/300)^(0.1) * exp(-(-4.5/T)),     # Found Match for Glover 38 CH + O --> CO + H and Umist 5310 CH + O --> CO + H with temp range [10, 2000]
         :k39 =>  n_H * 2.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 39 CH2 + H --> CH + H2 and Umist 5370 H + CH2 --> CH + H2 with temp range [10, 2500]
         :k40 =>  n_H * 1.33e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 40 CH2 + O --> CO + H + H and Umist 5232 CH2 + O --> CO + H + H with temp range [10, 2500]
         :k41 =>  n_H * 8e-11 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 41 CH2 + O --> CO + H2 and Umist 5231 CH2 + O --> CO + H2 with temp range [1900, 2600]
         :k42 =>  n_H * 2e-10 * (T/300)^(-0.12) * exp(-(0.0/T)),     # Found Match for Glover 42 C2 + O --> CO + C and Umist 5555 O + C2 --> CO + C with temp range [10, 8000]
         :k43 =>  n_H * 3.14e-13 * (T/300)^(2.7) * exp(-(3150.0/T)),     # Found Match for Glover 43 O + H2 --> OH + H and Umist 5361 H2 + O --> OH + H with temp range [297, 3532]
         :k44 =>  n_H * 6.99e-14 * (T/300)^(2.8) * exp(-(1950.0/T)),     # Found Match for Glover 44 OH + H --> O + H2 and Umist 5412 H + OH --> O + H2 with temp range [300, 2500]
         :k45 =>  n_H * 2.05e-12 * (T/300)^(1.52) * exp(-(1736.0/T)),     # Found Match for Glover 45 OH + H2 --> H2O + H and Umist 5362 H2 + OH --> H2O + H with temp range [250, 2581]
         :k46 =>  n_H * 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 46 OH + C --> CO + H and Umist 5199 C + OH --> CO + H with temp range [10, 300]
         :k47 =>  n_H * 3.69e-11 * (T/300)^(-0.27) * exp(-(12.9/T)),     # Found Match for Glover 47 OH + O --> O2 + H and Umist 5643 O + OH --> O2 + H with temp range [10, 500]
         #:k48 => FIX 48
         :k49 =>  n_H * 1.59e-11 * (T/300)^(1.2) * exp(-(9610.0/T)),     # Found Match for Glover 49 H2O + H --> H2 + OH and Umist 5382 H + H2O --> OH + H2 with temp range [250, 3000]
         :k50 =>  n_H * 2.61e-10 * (T/300)^(0.0) * exp(-(8156.0/T)),     # Found Match for Glover 50 O2 + H --> OH + O and Umist 5404 H + O2 --> OH + O with temp range [250, 4000]
         :k51 =>  n_H * 3.16e-10 * (T/300)^(0.0) * exp(-(21890.0/T)),     # Found Match for Glover 51 O2 + H2 --> OH + OH and Umist 5359 H2 + O2 --> OH + OH with temp range [300, 2500]
         :k52 =>  n_H * 5.56e-11 * (T/300)^(0.41) * exp(-(-26.9/T)),     # Found Match for Glover 52 O2 + C --> CO + O and Umist 5196 C + O2 --> CO + O with temp range [10, 8000]
         :k53 =>  n_H * 1.1e-10 * (T/300)^(0.5) * exp(-(77700.0/T)),     # Found Match for Glover 53 CO + H --> C + OH and Umist 5377 H + CO --> OH + C with temp range [2590, 41000]
         :k54 =>  n_H * 2.08e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 54 H2⁺ + H2 --> H3⁺ + H and Umist 2604 H2⁺ + H2 --> H3⁺ + H with temp range [10, 41000]
         # :k55 =>  Glover 55 H3⁺ + H --> H2⁺ + H2 had 0 product matches
         :k56 =>  n_H * 2.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 56 C + H2⁺ --> CH⁺ + H and Umist 2594 H2⁺ + C --> CH⁺ + H with temp range [10, 41000]
         :k57 =>  n_H * 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 57 C + H3⁺ --> CH⁺ + H2 and Umist 2852 H3⁺ + C --> CH⁺ + H2 with temp range [10, 41000]
         :k58 =>  n_H * 1e-10 * (T/300)^(0.0) * exp(-(4640.0/T)),     # Found Match for Glover 58 C⁺ + H2 --> CH⁺ + H and Umist 2618 H2 + C⁺ --> CH⁺ + H with temp range [154, 3000]
         :k59 =>  n_H * 9.06e-10 * (T/300)^(-0.37) * exp(-(29.1/T)),     # Found Match for Glover 59 CH⁺ + H --> C⁺ + H2 and Umist 3055 H + CH⁺ --> C⁺ + H2 with temp range [10, 1000]
         :k60 =>  n_H * 1.2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 60 CH⁺ + H2 --> CH2⁺ + H and Umist 2647 H2 + CH⁺ --> CH2⁺ + H with temp range [10, 41000]
         :k61 =>  n_H * 3.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 61 CH⁺ + O --> CO⁺ + H and Umist 2189 CH⁺ + O --> CO⁺ + H with temp range [10, 41000]
         :k62 =>  n_H * 1.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 62 CH2 + H⁺--> CH⁺ + H2 and Umist 2534 H⁺ + CH2 --> CH⁺ + H2 with temp range [10, 41000]
         :k63 =>  n_H * 1e-09 * (T/300)^(0.0) * exp(-(7080.0/T)),     # Found Match for Glover 63 CH2⁺ + H --> CH⁺ + H2 and Umist 3056 H + CH2⁺ --> CH⁺ + H2 with temp range [236, 300]
         :k64 =>  n_H * 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 64 CH2⁺ + H2 --> CH3⁺ + H and Umist 2648 H2 + CH2⁺ --> CH3⁺ + H with temp range [10, 41000]
         :k65 =>  n_H * 7.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 65 CH2⁺ + O --> HCO⁺ + H and Umist 2211 CH2⁺ + O --> HCO⁺ + H with temp range [10, 41000]
         :k66 =>  n_H * 7e-10 * (T/300)^(0.0) * exp(-(10560.0/T)),     # Found Match for Glover 66 CH3⁺ + H --> CH2⁺ + H2 and Umist 3057 H + CH3⁺ --> CH2⁺ + H2 with temp range [352, 41000]
         :k67 =>  n_H * 4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 67 CH3⁺ + O --> HCO⁺ + H2 and Umist 2311 CH3⁺ + O --> HCO⁺ + H2 with temp range [10, 41000]
         :k68 =>  n_H * 4.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 68 C2 + O⁺ --> CO⁺ + C and Umist 3851 O⁺ + C2 --> CO⁺ + C with temp range [10, 41000]
         :k69 =>  n_H * 1.7e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 69 O⁺ + H2 --> OH⁺ + H and Umist 2686 H2 + O⁺ --> OH⁺ + H with temp range [10, 41000]
         :k70 =>  n_H * 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 70 O + H2⁺ --> OH⁺ + H and Umist 2616 H2⁺ + O --> OH⁺ + H with temp range [10, 41000]
         :k71 =>  n_H * 7.98e-10 * (T/300)^(-0.16) * exp(-(1.4/T)),     # Found Match for Glover 71 O + H3⁺ --> OH⁺ + H2 and Umist 2951 H3⁺ + O --> OH⁺ + H2 with temp range [5, 400]
         :k72 =>  n_H * 1.3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 72 OH + H3⁺ --> H2O⁺ + H2 and Umist 2953 H3⁺ + OH --> H2O⁺ + H2 with temp range [10, 41000]
         :k73 =>  n_H * 7.7e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 73 OH + C⁺ --> CO⁺ + H and Umist 1646 C⁺ + OH --> CO⁺ + H with temp range [10, 41000]
         :k74 =>  n_H * 1.01e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 74 OH⁺ + H2 --> H2O⁺ + H and Umist 2689 H2 + OH⁺ --> H2O⁺ + H with temp range [10, 41000]
         :k75 =>  n_H * 6.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 75 H2O⁺ + H2 --> H3O⁺ + H and Umist 2662 H2 + H2O⁺ --> H3O⁺ + H with temp range [10, 41000]
         :k76 =>  n_H * 5.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 76 H2O + H3⁺ --> H3O⁺ + H2 and Umist 2904 H3⁺ + H2O --> H3O⁺ + H2 with temp range [10, 41000]
         :k77 =>  n_H * 9e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 77 H2O + C⁺ --> HCO⁺ + H and Umist 1613 C⁺ + H2O --> HCO⁺ + H with temp range [10, 41000]
         :k78 =>  n_H * 2.09e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 78 H2O + C⁺ --> HOC⁺ + H and Umist 1614 C⁺ + H2O --> HOC⁺ + H with temp range [10, 41000]
         :k79 =>  n_H * 1e-11 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 79 H3O⁺ + C --> HCO⁺ + H2 and Umist 2121 C + H3O⁺ --> HCO⁺ + H2 with temp range [10, 41000]
         :k80 =>  n_H * 3.42e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 80 O2 + C⁺ --> CO⁺ + O and Umist 1642 C⁺ + O2 --> CO⁺ + O with temp range [10, 41000]
         :k81 =>  n_H * 4.54e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 81 O2 + C⁺ --> CO + O⁺ and Umist 1643 C⁺ + O2 --> CO + O⁺ with temp range [10, 41000]
         :k82 =>  n_H * 9.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 82 O2 + CH2⁺ --> HCO⁺ + OH and Umist 2210 CH2⁺ + O2 --> HCO⁺ + OH with temp range [10, 41000]
         :k83 =>  n_H * 5.2e-11 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 83 O2⁺ + C --> CO⁺ + O and Umist 2132 C + O2⁺ --> CO⁺ + O with temp range [10, 41000]
         :k84 =>  n_H * 8.49e-10 * (T/300)^(0.07) * exp(-(5.2/T)),     # Found Match for Glover 84 CO + H3⁺ --> HOC⁺ + H2 and Umist 2896 H3⁺ + CO --> HOC⁺ + H2 with temp range [10, 400]
         :k85 =>  n_H * 1.36e-09 * (T/300)^(-0.14) * exp(-(-3.4/T)),     # Found Match for Glover 85 CO + H3⁺ --> HCO⁺ + H2 and Umist 2895 H3⁺ + CO --> HCO⁺ + H2 with temp range [10, 400]
         :k86 =>  n_H * 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 86 HCO⁺ + C --> CO + CH⁺ and Umist 2124 C + HCO⁺ --> CO + CH⁺ with temp range [10, 41000]
         :k87 =>  n_H * 2.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 87 HCO⁺ + H2O --> CO + H3O⁺ and Umist 2767 H2O + HCO⁺ --> CO + H3O⁺ with temp range [10, 41000]
         :k88 =>  n_H * 7.2e-15 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 88 H2 + He⁺ --> He + H2⁺ and Umist 458 H2 + He⁺ --> He + H2⁺ with temp range [10, 300]
         :k89 =>  n_H * 3.7e-14 * (T/300)^(0.0) * exp(-(35.0/T)),     # Found Match for Glover 89 H2 + He⁺ --> He + H + H⁺ and Umist 2677 H2 + He⁺ --> He + H⁺ + H with temp range [10, 300]
         :k90 =>  n_H * 1.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 90 CH + H⁺ --> CH⁺ + H and Umist 372 H⁺ + CH --> CH⁺ + H with temp range [10, 41000]
         :k91 =>  n_H * 1.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 91 CH2 + H⁺ --> CH2⁺ + H and Umist 359 H⁺ + CH2 --> CH2⁺ + H with temp range [10, 41000]
         :k92 =>  n_H * 7.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 92 CH2 + He⁺ --> C⁺ + He + H2 and Umist 3380 He⁺ + CH2 --> C⁺ + He + H2 with temp range [10, 41000]
         :k93 =>  n_H * 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 93 C2 + He⁺ --> C⁺ + C + He and Umist 3293 He⁺ + C2 --> C⁺ + C + He with temp range [10, 41000]
         :k94 =>  n_H * 2.1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 94 OH + H⁺ --> OH⁺ + H and Umist 408 H⁺ + OH --> OH⁺ + H with temp range [10, 41000]
         :k95 =>  n_H * 1.1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 95 OH + He⁺ --> O⁺ + He + H and Umist 3517 He⁺ + OH --> O⁺ + He + H with temp range [10, 41000]
         :k96 =>  n_H * 6.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 96 H2O + H⁺ --> H2O⁺ + H and Umist 378 H⁺ + H2O --> H2O⁺ + H with temp range [10, 300]
         :k97 =>  n_H * 2.04e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 97 H2O + He⁺ --> OH + He + H⁺ and Umist 3446 He⁺ + H2O --> OH + He + H⁺ with temp range [10, 41000]
         :k98 =>  n_H * 2.86e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 98 H2O + He⁺ --> OH⁺ + He + H and Umist 3445 He⁺ + H2O --> OH⁺ + He + H with temp range [10, 41000]
         :k99 =>  n_H * 6.05e-11 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 99 H2O + He⁺ --> H2O⁺ + He and Umist 518 He⁺ + H2O --> H2O⁺ + He with temp range [10, 41000]
         :k100 =>  n_H * 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 100 O2 + H⁺ --> O2⁺ + H and Umist 405 H⁺ + O2 --> O2⁺ + H with temp range [10, 300]
         :k101 =>  n_H * 3.3e-11 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 101 O2 + He⁺ --> O2⁺ + He and Umist 525 He⁺ + O2 --> O2⁺ + He with temp range [10, 41000]
         :k102 =>  n_H * 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 102 O2 + He⁺ --> O⁺ + O + He and Umist 3510 He⁺ + O2 --> O⁺ + O + He with temp range [10, 41000]
         :k103 =>  n_H * 5.2e-11 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 103 O2⁺ + C --> O2 + C⁺ and Umist 248 C + O2⁺ --> O2 + C⁺ with temp range [10, 41000]
         :k104 =>  n_H * .6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 104 CO + He⁺ --> C⁺ + O + He and Umist 3431 He⁺ + CO --> O + C⁺ + He with temp range [10, 41000]
         #:k105 => FIX
         :k106 =>  n_H * 7.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 106 CO⁺ + H --> CO + H⁺ and Umist 488 H + CO⁺ --> CO + H⁺ with temp range [10, 41000]
         :k107 =>  n_H * 7.51e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 107 C⁻ + H⁺ --> C + H and Umist 4116 C⁻ + H⁺ --> C + H with temp range [10, 41000]
         :k108 =>  n_H * 7.51e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 108 O⁻ + H⁺ --> O + H and Umist 4957 O⁻ + H⁺ --> O + H with temp range [10, 41000]
         :k109 =>  n_H * 7.51e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 109 He⁺ + H⁻ --> He + H and Umist 4932 H⁻ + He⁺ --> H + He with temp range [10, 41000]
         :k110 =>  n_H * 2.34e-08 * (T/300)^(-0.52) * exp(-(0.0/T)),     # Found Match for Glover 110 H3⁺ + e --> H2 + H and Umist 1277 H3⁺ + e --> H2 + H with temp range [10, 1000]
         :k111 =>  n_H * 4.36e-08 * (T/300)^(-0.52) * exp(-(0.0/T)),     # Found Match for Glover 111 H3⁺ + e --> H + H + H and Umist 1278 H3⁺ + e --> H + H + H with temp range [10, 1000]
         :k112 =>  n_H * 1.5e-07 * (T/300)^(-0.42) * exp(-(0.0/T)),     # Found Match for Glover 112 CH⁺ + e --> C + H and Umist 1158 CH⁺ + e --> C + H with temp range [10, 300]
         :k113 =>  n_H * 1.6e-07 * (T/300)^(-0.6) * exp(-(0.0/T)),     # Found Match for Glover 113 CH2⁺ + e --> CH + H and Umist 1161 CH2⁺ + e --> CH + H with temp range [10, 1000]
         :k114 =>  n_H * 4.03e-07 * (T/300)^(-0.6) * exp(-(0.0/T)),     # Found Match for Glover 114 CH2⁺ + e --> C + H + H and Umist 1160 CH2⁺ + e --> C + H + H with temp range [10, 1000]
         :k115 =>  n_H * 7.68e-08 * (T/300)^(-0.6) * exp(-(0.0/T)),     # Found Match for Glover 115 CH2⁺ + e --> C + H2 and Umist 1159 CH2⁺ + e --> C + H2 with temp range [10, 1000]
         :k116 =>  n_H * 7.75e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 116 CH3⁺ + e --> CH2 + H and Umist 1175 CH3⁺ + e --> CH2 + H with temp range [10, 300]
         :k117 =>  n_H * 1.95e-07 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 117 CH3⁺ + e --> CH + H2 and Umist 1176 CH3⁺ + e --> CH + H2 with temp range [10, 300]
         :k118 =>  n_H * 2e-07 * (T/300)^(-0.4) * exp(-(0.0/T)),     # Found Match for Glover 118 CH3⁺ + e --> CH + H + H and Umist 1177 CH3⁺ + e --> CH + H + H with temp range [10, 1000]
         :k119 =>  n_H * 3.75e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 119 OH⁺ + e --> O + H and Umist 1427 OH⁺ + e --> O + H with temp range [10, 300]
         :k120 =>  n_H * 3.05e-07 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 120 H2O⁺ + e --> O + H + H and Umist 1265 H2O⁺ + e --> O + H + H with temp range [10, 1000]
         :k121 =>  n_H * 3.9e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 121 H2O⁺ + e --> O + H2 and Umist 1264 H2O⁺ + e --> O + H2 with temp range [10, 1000]
         :k122 =>  n_H * 8.6e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 122 H2O⁺ + e --> OH + H and Umist 1266 H2O⁺ + e --> OH + H with temp range [10, 1000]
         :k123 =>  n_H * 7.09e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 123 H3O⁺ + e --> H + H2O and Umist 1293 H3O⁺ + e --> H2O + H with temp range [10, 1000]
         :k124 =>  n_H * 5.37e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 124 H3O⁺ + e --> OH + H2 and Umist 1295 H3O⁺ + e --> OH + H2 with temp range [10, 1000]
         :k125 =>  n_H * 3.05e-07 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 125 H3O⁺ + e --> OH + H + H and Umist 1296 H3O⁺ + e --> OH + H + H with temp range [10, 1000]
         :k126 =>  n_H * 5.6e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),     # Found Match for Glover 126 H3O⁺ + e --> O + H + H2 and Umist 1294 H3O⁺ + e --> O + H2 + H with temp range [10, 1000]
         :k127 =>  n_H * 1.95e-07 * (T/300)^(-0.7) * exp(-(0.0/T)),     # Found Match for Glover 127 O2⁺ + e --> O + O and Umist 1421 O2⁺ + e --> O + O with temp range [10, 300]
         :k128 =>  n_H * 2e-07 * (T/300)^(-0.48) * exp(-(0.0/T)),     # Found Match for Glover 128 CO⁺ + e --> C + O and Umist 1232 CO⁺ + e --> O + C with temp range [10, 41000]
         :k129 =>  n_H * 2.4e-07 * (T/300)^(-0.69) * exp(-(0.0/T)),     # Found Match for Glover 129 HCO⁺ + e --> CO + H and Umist 1349 HCO⁺ + e --> CO + H with temp range [10, 300]
         #:k130 => FIX
         :k131 =>  n_H * 1.1e-07 * (T/300)^(-1.0) * exp(-(0.0/T)),     # Found Match for Glover 131 HOC⁺ + e --> CO + H and Umist 1377 HOC⁺ + e --> CO + H with temp range [10, 300]
         :k132 =>  n_H * 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 132 H⁻ + C --> CH + e and Umist 69 H⁻ + C --> CH + e with temp range [10, 41000]
         :k133 =>  n_H * 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 133 H⁻ + O --> OH + e and Umist 80 H⁻ + O --> OH + e with temp range [10, 41000]
         :k134 =>  n_H * 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 134 H⁻ + OH --> H2O + e and Umist 81 H⁻ + OH --> H2O + e with temp range [10, 41000]
         :k135 =>  n_H * 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 135 C⁻ + H --> CH + e and Umist 84 H + C⁻ --> CH + e with temp range [10, 41000]
         :k136 =>  n_H * 1e-13 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 136 C⁻ + H2 --> CH2 + e and Umist 82 H2 + C⁻ --> CH2 + e with temp range [10, 41000]
         :k137 =>  n_H * 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 137 C⁻ + O --> CO + e and Umist 9 C⁻ + O --> CO + e with temp range [10, 41000]
         :k138 =>  n_H * 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 138 O⁻ + H --> OH + e and Umist 107 H + O⁻ --> OH + e with temp range [10, 41000]
         :k139 =>  n_H * 7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 139 O⁻ + H2 --> H2O + e and Umist 83 H2 + O⁻ --> H2O + e with temp range [10, 41000]
         :k140 =>  n_H * 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 140 O⁻ + C --> CO + e and Umist 57 C + O⁻ --> CO + e with temp range [10, 41000]
         # :k141 =>  Glover 141 H2 + H⁺ --> H3⁺ had 0 product matches
         :k142 =>  n_H * 2.25e-15 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 142 C + e --> C⁻  and Umist 6154 C + e --> C⁻ with temp range [10, 1000]
         :k143 =>  n_H * 1e-17 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 143 C + H --> CH  and Umist 6119 H + C --> CH with temp range [10, 300]
         :k144 =>  n_H * 1e-17 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 144 C + H2 --> CH2  and Umist 6101 H2 + C --> CH2 with temp range [10, 300]
         :k145 =>  n_H * 4.36e-18 * (T/300)^(0.35) * exp(-(161.3/T)),     # Found Match for Glover 145 C + C --> C2  and Umist 6070 C + C --> C2 with temp range [10, 41000]
         :k146 =>  n_H * 4.69e-19 * (T/300)^(1.52) * exp(-(-50.5/T)),     # Found Match for Glover 146 C + O --> CO  and Umist 6073 C + O --> CO with temp range [10, 14700]
         :k147 =>  n_H * 1.7e-17 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 147 C⁺ + H --> CH⁺  and Umist 6116 H + C⁺ --> CH⁺ with temp range [10, 300]
         :k148 =>  n_H * 2e-16 * (T/300)^(-1.3) * exp(-(23.0/T)),     # Found Match for Glover 148 C⁺ + H2 --> CH2⁺  and Umist 6095 H2 + C⁺ --> CH2⁺ with temp range [10, 300]
         :k149 =>  n_H * 3.14e-18 * (T/300)^(-0.15) * exp(-(68.0/T)),     # Found Match for Glover 149 C⁺ + O --> CO⁺  and Umist 6052 C⁺ + O --> CO⁺ with temp range [10, 13900]
         :k150 =>  n_H * 1.5e-15 * (T/300)^(0.0) * exp(-(0.0/T)),     # Found Match for Glover 150 O + e --> O⁻  and Umist 6156 O + e --> O⁻ with temp range [10, 41000]
         :k151 =>  n_H * 9.9e-19 * (T/300)^(-0.38) * exp(-(0.0/T)),     # Found Match for Glover 151 O + H --> OH  and Umist 6120 H + O --> OH with temp range [10, 300]
         :k152 =>  n_H * 4.9e-20 * (T/300)^(1.58) * exp(-(0.0/T)),     # Found Match for Glover 152 O + O --> O2  and Umist 6129 O + O --> O2 with temp range [10, 300]
         :k153 =>  n_H * 5.26e-18 * (T/300)^(-5.22) * exp(-(90.0/T)),     # Found Match for Glover 153 OH + H --> H2O  and Umist 6121 H + OH --> H2O with temp range [20, 300]
         
         
         :k154 => n_H^2 * 1.32e-32 * (T/300)^-0.38,   # RATE FROM Glover reaction 154 H + H + H --> H2 + H
         :k155 => n_H^2 * 2.8e-31 * (T)^-0.6,         # RATE FROM Glover reaction 155 H + H + H2 --> H2 + H2
         :k156 => n_H^2 * 6.9e-32 * (T)^(-0.4),       # RATE FROM Glover reaction 156 H + H + He --> H2 + He
         :k157 => n_H^2 * 5.99e-33 * (T/5000)^(-1.6), # RATE FROM Glover reaction 157 C + C + M --> C2 + M
         :k158 => n_H^2 * 6.16e-29 * (T/300)^(-3.08), # RATE FROM Glover reaction 158 C + O + M --> CO + M
         :k159 => n_H^2 * 100 * cr_ion_rate * 960,    # RATE FROM Glover reaction 159 C⁺ + O + M --> CO⁺ + M
         :k160 => n_H^2 * 100 * cr_ion_rate * 960,    # RATE FROM Glover reaction 160 C + O⁺ + M --> CO⁺ + M
         :k161 => n_H^2 * 4.33e-32 * (T/300)^(-1),    # RATE FROM Glover reaction 161 O + H + M --> OH + M
         :k162 => n_H^2 * 2.56e-31 * (T/300)^(-2),    # RATE FROM Glover reaction 162 OH + H + M --> H2O + M
         :k163 => n_H^2 * 9.2e-34 * (T/300)^(-1),     # RATE FROM Glover reaction 163 O + O + M --> O2 + M
         

         :k164 =>  n_H * 1.09e-11 * (T/300)^(-2.19) * exp(-(165.1/T)),     # Found Match for Glover 164 O + CH --> HCO⁺ + e and Umist 64 CH + O --> HCO⁺ + e with temp range [10, 2500]
         # :k165 =>  Glover 165 H + H --> H2 had 0 product matches
         :k166 =>  1.43e-07 * exp(-(0.5*Av)),     # Found Match for Glover 166 H⁻ --> H + e and Umist 5902 (PH) H⁻ --> H + e with temp range [10, 41000]
         :k167 =>  5.7e-10 * exp(-(2.4*Av)),     # Found Match for Glover 167 H2⁺ --> H + H⁺ and Umist 5903 H2⁺ --> H⁺ + H with temp range [10, 41000]
         :k168 => 1.3e-18,     # Found Match for Glover 168 H2 --> H + H and Umist 732 H2 --> H + H with temp range [10, 41000]
         :k169 =>  5e-15 * exp(-(1.8*Av)),     # Found Match for Glover 169 H3⁺ --> H2 + H⁺ and Umist 5923 H3⁺ --> H2 + H⁺ with temp range [10, 41000]
         :k170 =>  5e-15 * exp(-(2.3*Av)),     # Found Match for Glover 170 H3⁺ --> H2⁺ + H and Umist 5922 H3⁺ --> H2⁺ + H with temp range [10, 41000]
         :k171 =>  3.1e-10 * exp(-(3.3*Av)),     # Found Match for Glover 171 C --> C⁺ + e and Umist 5827 C --> C⁺ + e with temp range [10, 41000]
         :k172 =>  4.9e-08 * exp(-(0.5*Av)),     # Found Match for Glover 172 C⁻ --> C + e and Umist 5706 C⁻ --> C + e with temp range [10, 41000]
         :k173 =>  9.2e-10 * exp(-(1.7*Av)),     # Found Match for Glover 173 CH --> C + H and Umist 5887 CH --> C + H with temp range [10, 41000]
         :k174 =>  7.6e-10 * exp(-(3.3*Av)),     # Found Match for Glover 174 CH --> CH⁺ + e and Umist 5888 CH --> CH⁺ + e with temp range [10, 41000]
         :k175 =>  3.3e-10 * exp(-(2.9*Av)),     # Found Match for Glover 175 CH⁺ --> C + H⁺ and Umist 5831 CH⁺ --> C + H⁺ with temp range [10, 41000]
         :k176 =>  5.8e-10 * exp(-(2.0*Av)),     # Found Match for Glover 176 CH2 --> CH + H and Umist 5837 CH2 --> CH + H with temp range [10, 41000]
         :k177 =>  1e-09 * exp(-(2.3*Av)),     # Found Match for Glover 177 CH2 --> CH2⁺ + e and Umist 5836 CH2 --> CH2⁺ + e with temp range [10, 41000]
         :k178 =>  4.67e-11 * exp(-(2.2*Av)),     # Found Match for Glover 178 CH2⁺ --> CH⁺ + H and Umist 5834 CH2⁺ --> CH⁺ + H with temp range [10, 41000]
         :k179 =>  1e-09 * exp(-(1.7*Av)),     # Found Match for Glover 179 CH3⁺ --> CH2⁺ + H and Umist 5854 CH3⁺ --> CH2⁺ + H with temp range [10, 41000]
         :k180 =>  1e-09 * exp(-(1.7*Av)),     # Found Match for Glover 180 CH3⁺ --> CH⁺ + H2 and Umist 5853 CH3⁺ --> CH⁺ + H2 with temp range [10, 41000]
         :k181 =>  2.4e-10 * exp(-(2.6*Av)),     # Found Match for Glover 181 C2 --> C + C and Umist 5736 C2 --> C + C with temp range [10, 41000]
         :k182 =>  1.09e-08 * exp(-(0.5*Av)),     # Found Match for Glover 182 O⁻ --> O + e and Umist 5985 O⁻ --> O + e with temp range [10, 41000]
         :k183 =>  3.9e-10 * exp(-(2.2*Av)),     # Found Match for Glover 183 OH --> O + H and Umist 5997 OH --> O + H with temp range [10, 41000]
         :k184 =>  1.6e-12 * exp(-(3.1*Av)),     # Found Match for Glover 184 OH --> OH⁺ + e and Umist 5998 OH --> OH⁺ + e with temp range [10, 41000]
         # :k185 => NO match for Glover reaction 185 OH⁺ --> O + H⁺
         :k186 =>  8e-10 * exp(-(2.2*Av)),     # Found Match for Glover 186 H2O --> OH + H and Umist 5915 H2O --> OH + H with temp range [10, 41000]
         :k187 =>  3.1e-11 * exp(-(3.9*Av)),     # Found Match for Glover 187 H2O --> H2O⁺ + e and Umist 5914 H2O --> H2O⁺ + e with temp range [10, 41000]
         # :k188 => NO match for Glover reaction 188 H2O⁺ --> H2⁺ + O
         # :k189 => NO match for Glover reaction 189 H2O⁺ --> H⁺ + OH
         # :k190 => NO match for Glover reaction 190 H2O⁺ --> O⁺ + H2
         :k191 =>  1e-12 * exp(-(2.0*Av)),     # Found Match for Glover 191 H2O⁺ --> OH⁺ + H and Umist 5912 H2O⁺ --> OH⁺ + H with temp range [10, 41000]
         # :k192 => NO match for Glover reaction 192 H3O⁺ --> H⁺ + H2O
         # :k193 => NO match for Glover reaction 193 H3O⁺ --> H2⁺ + OH
         # :k194 => NO match for Glover reaction 194 H3O⁺ --> H2O⁺ + H
         # :k195 => NO match for Glover reaction 195 H3O⁺ --> OH⁺ + H2
         :k196 =>  7.6e-11 * exp(-(3.9*Av)),     # Found Match for Glover 196 O2 --> O2⁺ + e and Umist 5988 O2 --> O2⁺ + e with temp range [10, 41000]
         :k197 =>  7.9e-10 * exp(-(2.1*Av)),     # Found Match for Glover 197 O2 --> O + O and Umist 5989 O2 --> O + O with temp range [10, 41000]
         :k198 =>  2e-10 * exp(-(3.5*Av)),     # Found Match for Glover 198 CO --> C + O and Umist 5894 CO --> O + C with temp range [10, 41000]
         :k199 => 5.98e-18,     # Found Match for Glover 199 H --> H⁺ + e and Umist 733 (CP) H --> H⁺ + e with temp range [10, 41000]
         :k200 => 6.5e-18,     # Found Match for Glover 200 He --> He⁺ + e and Umist 734 (CP) He --> He⁺ + e with temp range [10, 41000]
         :k201 => 2.86e-19,     # Found Match for Glover 201 H2 --> H⁺ + H + e and Umist 730 H2 --> H⁺ + H + e with temp range [10, 41000]
         :k202 => 1.3e-18,     # Found Match for Glover 202 H2 --> H + H and Umist 732 H2 --> H + H with temp range [10, 41000]
         :k203 => 3.9e-21,     # Found Match for Glover 203 H2 --> H⁺ + H⁻ and Umist 729 H2 --> H⁺ + H⁻ with temp range [10, 41000]
         :k204 => 1.2e-17,     # Found Match for Glover 204 H2 --> H2⁺ + e and Umist 731 H2 --> H2⁺ + e with temp range [10, 41000]
         :k205 => 2.3e-17,     # Found Match for Glover 205 C --> C⁺ + e and Umist 726 (CP) C --> C⁺ + e with temp range [10, 41000]
         :k206 => 3.4e-17,     # Found Match for Glover 206 O --> O⁺ + e and Umist 736 (CP) O --> O⁺ + e with temp range [10, 41000]
         :k207 => 3.9e-17,     # Found Match for Glover 207 CO --> CO⁺ + e and Umist 727 CO --> CO⁺ + e with temp range [10, 41000]
         :k208 =>  1.3e-17 * (T/300)^(0.0) * (255.0/(1-omega)),     # Found Match for Glover 208 C --> C⁺ + e and Umist 823 (CR) C --> C⁺ + e with temp range [10, 41000]
         :k209 =>  1.3e-17 * (T/300)^(0.0) * (365.0/(1-omega)),     # Found Match for Glover 209 CH --> C + H and Umist 871 CH --> C + H with temp range [10, 41000]
         :k210 =>  1.3e-17 * (T/300)^(0.0) * (88.0/(1-omega)),     # Found Match for Glover 210 CH⁺ --> C⁺ + H and Umist 827 CH⁺ --> C⁺ + H with temp range [10, 41000]
         :k211 =>  1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),     # Found Match for Glover 211 CH2 --> CH2⁺ + e and Umist 829 CH2 --> CH2⁺ + e with temp range [10, 41000]
         :k212 =>  1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),     # Found Match for Glover 212 CH2 --> CH + H and Umist 830 CH2 --> CH + H with temp range [10, 41000]
         :k213 =>  1.3e-17 * (T/300)^(0.0) * (119.5/(1-omega)),     # Found Match for Glover 213 C2 --> C + C and Umist 745 C2 --> C + C with temp range [10, 41000]
         :k214 =>  1.3e-17 * (T/300)^(0.0) * (254.5/(1-omega)),     # Found Match for Glover 214 OH --> O + H and Umist 955 OH --> O + H with temp range [10, 41000]
         :k215 =>  1.3e-17 * (T/300)^(0.0) * (485.5/(1-omega)),     # Found Match for Glover 215 H2O --> OH + H and Umist 889 H2O --> OH + H with temp range [10, 41000]
         :k216 =>  1.3e-17 * (T/300)^(0.0) * (375.5/(1-omega)),     # Found Match for Glover 216 O2 --> O + O and Umist 948 O2 --> O + O with temp range [10, 41000]
         :k217 =>  1.3e-17 * (T/300)^(0.0) * (58.5/(1-omega)),     # Found Match for Glover 217 O2 --> O2⁺ + e and Umist 947 O2 --> O2⁺ + e with temp range [10, 41000]
         :k218 =>  1.3e-17 * (T/300)^(1.17) * (105.0/(1-omega)),     # Found Match for Glover 218 CO --> C + O and Umist 876 CO --> O + C with temp range [10, 41000]
         )



### Tempurature Range Flags ###
T_low_bounds = [10  10 200 10 10 10 -1 3400 1833 2803 -1 10 -1 -1 -1 -1 10 10 -1 10  10 -1 -1 10 10 -1 -1 -1 10 -1 1696 20 -1 300 300 300 10 10 10 10 1900 10 297 300 250 10 10 -1 250 250 300 10 2590 10 -1 10 10 154 10 10 10 10 236 10 10 352 10 10 10 10 5 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 -1 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 -1 10 10 10 10 10 10 10 10 10 10 -1 10 10 10 10 10 10 10 10 10 10 10 20 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 10 -1 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 -1 10 10 -1 -1 -1 10 -1 -1 -1 -1 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10]
                    
T_upp_bounds = [41000 100 32000 41000 41000 300 1e10 41000 41000 41000 1e10 20000 1e10 1e10 1e10 1e10 1000 41000 1e10 41000 41000 1e10 1e10 41000 41000 1e10 1e10 1e10 300 1e10 41000 41000 1e10 2500 2000 2500 300 2000 2500 2500 2600 8000 3532 2500 2581 300 500 1e6 3000 4000 2500 8000 41000 41000 1e10 41000 41000 3000 1000 41000 41000 41000 300 41000 41000 41000 41000 41000 41000 41000 400 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 400 400 41000 41000 300 300 41000 41000 41000 41000 41000 41000 300 41000 41000 41000 300 41000 41000 41000 41000 1e6 41000 41000 41000 41000 1000 1000 300 1000 1000 1000 300 300 1000 300 1000 1000 1000 1000 1000 1000 1000 300 41000 300 1e6 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 1e10 1000 300 300 41000 14700 300 300 13900 41000 300 300 300 1e10 1e10 1e10 1e10 1e10 1e10 1e10 1e10 1e10 1e10 2500 1e10 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 1e10 41000 41000 1e10 1e10 1e10 41000 1e10 1e10 1e10 1e10 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000]
    
for i = 1:length(T_low_bounds)
    if T < T_low_bounds[i]
        print("\nTempurature below the tempurature range [", T_low_bounds[i], ", ", T_upp_bounds[i], "] for reaction number ", i ) 
    elseif T > T_upp_bounds[i]
        print("\nTempurature higher than the tempurature range [", T_low_bounds[i], ", ", T_upp_bounds[i], "] for reaction number ", i ) 
    end
end


# %% Network admin things
# FIX M
@variables t 
@species H⁻(t) H2⁺(t) H3⁺(t) CH⁺(t) CH2⁺(t) OH⁺(t) H3O⁺(t) CO⁺(t) HOC⁺(t) O⁻(t) C⁻(t) O2⁺(t) e(t) H⁺(t) H(t) H2(t) He(t) He⁺(t) C(t) C⁺(t) O(t) O⁺(t) OH(t) CO(t) C2(t) O2(t) H2O⁺(t) CH(t) CH2(t) CH3⁺(t) M(t)
@parameters T Av n_H cr_ion_rate omega k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18 k19 k20 k21 k22 k23 k24 k25 k26 k27 k28 k29 k30 k31 k32 k33 k34 k35 k36 k37 k38 k39 k40 k41 k42 k43 k44 k45 k46 k47 k48 k49 k50 k51 k52 k53 k54 k55 k56 k57 k58 k59 k60 k61 k62 k63 k64 k65 k66 k67 k68 k69 k70 k71 k72 k73 k74 k75 k76 k77 k78 k79 k80 k81 k82 k83 k84 k85 k86 k87 k88 k89 k90 k91 k92 k93 k94 k95 k96 k97 k98 k99 k100 k101 k102 k103 k104 k105 k106 k107 k108 k109 k110 k111 k112 k113 k114 k115 k116 k117 k118 k119 k120 k121 k122 k123 k124 k125 k126 k127 k128 k129 k130 k131 k132 k133 k134 k135 k136 k137 k138 k139 k140 k141 k142 k143 k144 k145 k146 k147 k148 k149 k150 k151 k152 k153 k154 k155 k156 k157 k158 k159 k160 k161 k162 k163 k164 k165 k166 k167 k168 k169 k170 k171 k172 k173 k174 k175 k176 k177 k178 k179 k180 k181 k182 k183 k184 k185 k186 k187 k188 k189 k190 k191 k192 k193 k194 k195 k196 k197 k198 k199 k200 k201 k202 k203 k204 k205 k206 k207 k208 k209 k210 k211 k212 k213 k214 k215 k216 k217 k218 


### Reaction Equations ###
reaction_equations = [
    (@reaction k1, H + e --> H⁻),  
    (@reaction k2, H⁻ + H --> H2 + e),
    (@reaction k3, H + H⁺ --> H2⁺ ), 
    (@reaction k4, H + H2⁺ --> H2 + H⁺), 
    (@reaction k5, H⁻ + H⁺--> H + H), 
    (@reaction k6, H2⁺ + e --> H + H), 
    #(@reaction k7, H2 + H⁺ --> H2⁺ + H), 
    (@reaction k8, H2 + e --> H + H + e), 
    (@reaction k9, H2 + H --> H + H + H), 
    (@reaction k10, H2 + H2 --> H2 + H + H), 
    #(@reaction k11, H + e --> H⁺ + e + e), 
    (@reaction k12, H⁺ + e --> H ),
    #(@reaction k13, H⁻ + e --> H + e + e), 
    #(@reaction k14, H⁻ + H --> H + H + e), 
    #(@reaction k15, H⁻ + H⁺--> H2⁺ + e),
    #(@reaction k16, He + e --> He⁺ + e + e), 
    (@reaction k17, He⁺ + e --> He ), 
    (@reaction k18, He⁺ + H --> He + H⁺), 
    #(@reaction k19, He + H⁺--> He⁺ + H), 
    
    (@reaction k20, C⁺ + e --> C ),
    (@reaction k20_backwards, C --> C⁺ + e ),

    (@reaction k21, O⁺ + e --> O ), 
    #(@reaction k22, C + e --> C⁺ + e + e), 
    #(@reaction k23, O + e --> O⁺ + e + e), 
    (@reaction k24, O⁺ + H --> O + H⁺), 
    (@reaction k25, O + H⁺--> O⁺ + H), 
    #(@reaction k26, O + He⁺ --> O⁺ + He), 
    #(@reaction k27, C + H⁺ --> C⁺ + H), 
    #(@reaction k28, C⁺ + H --> C + H⁺), 
    (@reaction k29, C + He⁺ --> C⁺ + He), 
    #(@reaction k30, H2 + He --> H + H + He), 
    (@reaction k31, OH + H --> O + H + H), 
    (@reaction k32, HOC⁺ + H2 --> HCO⁺ + H2), 
    #(@reaction k33, HOC⁺ + CO --> HCO⁺ + CO), 
    (@reaction k34, C + H2 --> CH + H), 
    (@reaction k35, CH + H --> C + H2), 
    (@reaction k36, CH + H2 --> CH2 + H), 
    (@reaction k37, CH + C --> C2 + H), 
    (@reaction k38, CH+O --> CO+H), # Nelson Match Neutral-Neutral reaction 1, Nelson has 2e-10
    (@reaction k39, CH2 + H --> CH + H2), 
    (@reaction k40, CH2 + O --> CO + H + H), 
    (@reaction k41, CH2 + O --> CO + H2), 
    (@reaction k42, C2 + O --> CO + C), 
    (@reaction k43, O + H2 --> OH + H), 
    (@reaction k44, OH + H --> O + H2), 
    (@reaction k45, OH + H2 --> H2O + H), 
    (@reaction k46, OH + C --> CO + H), # Nelson Match Neutral-Neutral reaction 2, Nelson has 5.8e-12 * T^(0.5)
    (@reaction k47, OH + O --> O2 + H), 
    #(@reaction k48, OH + OH --> H2O + H), 
    (@reaction k49, H2O + H --> H2 + OH), 
    (@reaction k50, O2 + H --> OH + O), 
    (@reaction k51, O2 + H2 --> OH + OH), 
    (@reaction k52, O2 + C --> CO + O), 
    (@reaction k53, CO + H --> C + OH), 
    (@reaction k54, H2⁺ + H2 --> H3⁺ + H), 
    #(@reaction k55, H3⁺ + H --> H2⁺ + H2), 
    (@reaction k56, C + H2⁺ --> CH⁺ + H), 
    (@reaction k57, C + H3⁺ --> CH⁺ + H2), # Nelson Match Ion-Molecule reaction 1
    (@reaction k58, C⁺ + H2 --> CH⁺ + H), # Nelson Match Ion-Molecule reaction 6, Nelson has 4e-16
    (@reaction k59, CH⁺ + H --> C⁺ + H2), 
    (@reaction k60, CH⁺ + H2 --> CH2⁺ + H), 
    (@reaction k61, CH⁺ + O --> CO⁺ + H), 
    (@reaction k62, CH2 + H⁺--> CH⁺ + H2), 
    (@reaction k63, CH2⁺ + H --> CH⁺ + H2), 
    (@reaction k64, CH2⁺ + H2 --> CH3⁺ + H), 
    (@reaction k65, CH2⁺ + O --> HCO⁺ + H), 
    (@reaction k66, CH3⁺ + H --> CH2⁺ + H2), 
    (@reaction k67, CH3⁺ + O --> HCO⁺ + H2), 
    (@reaction k68, C2 + O⁺ --> CO⁺ + C), 
    (@reaction k69, O⁺ + H2 --> OH⁺ + H), 
    (@reaction k70, O + H2⁺ --> OH⁺ + H), 
    (@reaction k71, O + H3⁺ --> OH⁺ + H2), # Nelson Match Ion-Molecule reaction 2 Nelson has 8e-10
    (@reaction k72, OH + H3⁺ --> H2O⁺ + H2), 
    (@reaction k73, OH + C⁺ --> CO⁺ + H), 
    (@reaction k74, OH⁺ + H2 --> H2O⁺ + H), 
    (@reaction k75, H2O⁺ + H2 --> H3O⁺ + H), 
    (@reaction k76, H2O + H3⁺ --> H3O⁺ + H2), 
    (@reaction k77, H2O + C⁺ --> HCO⁺ + H), 
    (@reaction k78, H2O + C⁺ --> HOC⁺ + H), 
    (@reaction k79, H3O⁺ + C --> HCO⁺ + H2), 
    (@reaction k80, O2 + C⁺ --> CO⁺ + O), 
    (@reaction k81, O2 + C⁺ --> CO + O⁺), 
    (@reaction k82, O2 + CH2⁺ --> HCO⁺ + OH), 
    (@reaction k83, O2⁺ + C --> CO⁺ + O), 
    (@reaction k84, CO + H3⁺ --> HOC⁺ + H2), 
    (@reaction k85, CO + H3⁺ --> HCO⁺ + H2), # Nelson Match Ion-Molecule reaction 3
    (@reaction k86, HCO⁺ + C --> CO + CH⁺), 
    (@reaction k87, HCO⁺ + H2O --> CO + H3O⁺), 
    (@reaction k88, H2 + He⁺  -->  He + H2⁺), # Nelson Match Ion-Molecule reaction 4, Nelson has 7e-15
    (@reaction k89, H2 + He⁺ --> He + H + H⁺), # Nelson Match Ion-Molecule reaction 4, Nelson has 7e-15
    (@reaction k90, CH + H⁺ -->  CH⁺ + H), 
    (@reaction k91, CH2 + H⁺ -->  CH2⁺ + H), 
    (@reaction k92, CH2 + He⁺ -->  C⁺ + He + H2), 
    (@reaction k93, C2 + He⁺ -->  C⁺ + C + He), 
    (@reaction k94, OH + H⁺ -->  OH⁺ + H), 
    (@reaction k95, OH + He⁺ -->  O⁺ + He + H), 
    (@reaction k96, H2O + H⁺ -->  H2O⁺ + H), 
    (@reaction k97, H2O + He⁺ -->  OH + He + H⁺), 
    (@reaction k98, H2O + He⁺ -->  OH⁺ + He + H), 
    (@reaction k99, H2O + He⁺ -->  H2O⁺ + He), 
    (@reaction k100, O2 + H⁺ -->  O2⁺ + H), 
    (@reaction k101, O2 + He⁺ -->  O2⁺ + He), 
    (@reaction k102, O2 + He⁺ -->  O⁺ + O + He), 
    (@reaction k103, O2⁺ + C  -->  O2 + C⁺), 
    (@reaction k104, CO + He⁺ -->  C⁺ + O + He), # Nelson Match Ion-Molecule reaction 5, Nelson has 1.6e-9
    #(@reaction k105, CO + He⁺ -->  C + O⁺ + He), 
    (@reaction k106, CO⁺ + H  -->  CO + H⁺), 
    (@reaction k107, C⁻ + H⁺ -->  C + H), 
    (@reaction k108, O⁻ + H⁺ -->  O + H), 
    (@reaction k109, He⁺ + H⁻  -->  He + H), 
    (@reaction k110, H3⁺ + e --> H2 + H), 
    (@reaction k111, H3⁺ + e --> H + H + H), 
    (@reaction k112, CH⁺ + e --> C + H), 
    (@reaction k113, CH2⁺ + e --> CH + H), 
    (@reaction k114, CH2⁺ + e --> C + H + H), 
    (@reaction k115, CH2⁺ + e --> C + H2), 
    (@reaction k116, CH3⁺ + e --> CH2 + H), 
    (@reaction k117, CH3⁺ + e --> CH + H2), 
    (@reaction k118, CH3⁺ + e --> CH + H + H), 
    (@reaction k119, OH⁺ + e --> O + H), 
    (@reaction k120, H2O⁺ + e --> O + H + H), 
    (@reaction k121, H2O⁺ + e --> O+H2), 
    (@reaction k122, H2O⁺ + e --> OH + H), 
    (@reaction k123, H3O⁺ + e --> H + H2O), 
    (@reaction k124, H3O⁺ + e --> OH + H2), 
    (@reaction k125, H3O⁺ + e --> OH + H + H), 
    (@reaction k126, H3O⁺ + e --> O + H + H2), 
    (@reaction k127, O2⁺ + e --> O + O), 
    (@reaction k128, CO⁺ + e --> C + O), 
    (@reaction k129, HCO⁺ + e --> CO + H), 
    #(@reaction k130, HCO⁺ + e --> OH + C), 
    (@reaction k131, HOC⁺ + e --> CO + H), 
    (@reaction k132, H⁻ + C --> CH + e), 
    (@reaction k133, H⁻ + O --> OH + e), 
    (@reaction k134, H⁻ + OH --> H2O + e), 
    (@reaction k135, C⁻ + H --> CH + e), 
    (@reaction k136, C⁻ + H2 --> CH2 + e), 
    (@reaction k137, C⁻ + O -->  CO + e), 
    (@reaction k138, O⁻ + H  -->  OH + e), 
    (@reaction k139, O⁻ + H2 -->  H2O + e), 
    (@reaction k140, O⁻ + C --> CO + e), 
    #(@reaction k141, H2 + H⁺ --> H3⁺), 
    (@reaction k142, C + e --> C⁻ ), 
    (@reaction k143, C + H --> CH ), 
    (@reaction k144, C + H2 --> CH2 ), 
    (@reaction k145, C + C --> C2 ), 
    (@reaction k146, C + O  --> CO ), 
    (@reaction k147, C⁺ + H --> CH⁺ ), 
    (@reaction k148, C⁺ + H2 --> CH2⁺ ), 
    (@reaction k149, C⁺ + O --> CO⁺ ), 
    (@reaction k150, O + e --> O⁻ ), 
    (@reaction k151, O+H --> OH ), 
    (@reaction k152, O + O --> O2 ), 
    (@reaction k153, OH + H --> H2O ), 
    (@reaction k154, H + H + H --> H2 + H), 
    (@reaction k155, H + H + H2 --> H2 + H2), 
    (@reaction k156, H + H + He --> H2 + He), 
    (@reaction k157, C + C + M --> C2 + M), 
    (@reaction k158, C + O + M --> CO + M), 
    (@reaction k159, C⁺ + O + M --> CO⁺ + M), # FIX: these rely on k210, which doesn't exist?
    (@reaction k160, C + O⁺ + M --> CO⁺ + M), # FIX: these rely on k210, which doesn't exist?
    (@reaction k161, O + H + M --> OH + M), 
    (@reaction k162, OH + H + M --> H2O + M), 
    (@reaction k163, O + O + M --> O2 + M), 
    (@reaction k164, O + CH -->  HCO⁺ + e), 
    #(@reaction k165, H+H --> H2), # FIX H is H(s) here?
    

    # Photochemical Reactions (R166) start below
    (@reaction k166, H⁻ --> H + e), 
    (@reaction k167, H2⁺ --> H + H⁺), 
    (@reaction k168, H2   -->  H+H), 
    (@reaction k169, H3⁺ --> H2 + H⁺),
    (@reaction k170, H3⁺ --> H2⁺ + H), 
    (@reaction k171, C -->  C⁺ + e), 
    (@reaction k172, C⁻ --> C + e), 
    (@reaction k173, CH --> C + H), 
    (@reaction k174, CH --> CH⁺ + e), 
    (@reaction k175, CH⁺ --> C + H⁺), 
    (@reaction k176, CH2 --> CH + H), 
    (@reaction k177, CH2 --> CH2⁺ + e), 
    (@reaction k178, CH2⁺ --> CH⁺ + H), 
    (@reaction k179, CH3⁺ --> CH2⁺ + H), 
    (@reaction k180, CH3⁺ --> CH⁺ + H2), 
    (@reaction k181, C2 --> C + C), 
    (@reaction k182, O⁻ --> O + e),
    (@reaction k183, OH --> O + H), 
    (@reaction k184, OH --> OH⁺ + e), 
    #(@reaction k185, OH⁺  --> O + H⁺), 
    (@reaction k186, H2O --> OH + H), 
    (@reaction k187, H2O --> H2O⁺ + e), 
    #(@reaction k188, H2O⁺ --> H2⁺ + O), 
    #(@reaction k189, H2O⁺ --> H⁺ + OH), 
    #(@reaction k190, H2O⁺ --> O⁺ + H2), 
    (@reaction k191, H2O⁺ --> OH⁺ + H), 
    #(@reaction k192, H3O⁺ --> H⁺ + H2O), 
    #(@reaction k193, H3O⁺ --> H2⁺ + OH), 
    #(@reaction k194, H3O⁺ --> H2O⁺ + H), 
    #(@reaction k195, H3O⁺ --> OH⁺ + H2), 
    (@reaction k196, O2 --> O2⁺ + e), 
    (@reaction k197, O2 --> O + O), 
    (@reaction k198, CO   --> C + O), #from UMIST


    # Cosmic Ray Reactions (R199) start below 
    (@reaction k199, H --> H⁺ + e), 
    (@reaction k200, He --> He⁺ + e), # cr_rate = 6e-18, Nelson match CR Ionization reaction 2
    (@reaction k201, H2 --> H⁺ + H + e), 
    (@reaction k202, H2 --> H + H), 
    (@reaction k203, H2 --> H⁺ + H⁻), 
    (@reaction k204, H2 --> H2⁺ + e),  # cr_rate = 6e-18, Nelson match CR Ionization reaction 1
    (@reaction k205, C --> C⁺ + e), 
    (@reaction k206, O --> O⁺ + e), 
    (@reaction k207, CO --> CO⁺ + e), 
    (@reaction k208, C --> C⁺ + e), 
    (@reaction k209, CH --> C + H), 
    (@reaction k210, CH⁺ --> C⁺ + H), 
    (@reaction k211, CH2 --> CH2⁺ + e), 
    (@reaction k212, CH2 --> CH + H), 
    (@reaction k213, C2 --> C + C), 
    (@reaction k214, OH --> O + H), 
    (@reaction k215, H2O --> OH + H), 
    (@reaction k216, O2 --> O + O), 
    (@reaction k217, O2 --> O2⁺ + e), 
    (@reaction k218, CO --> C + O) # this is a cosmic ray reaction, NOT photoreaction
]

print("\nCheckpoint 3: Finished reading reaction equations")
### Turn the Network into a system of ODEs ###
@named system = ReactionSystem(reaction_equations, t)
print("\nCheckpoint 4: Finished creating the reaction system")
sys = convert(ODESystem, complete(system))
print("\nCheckpoint 5: Finished converting to an ODE System")
ssys = structural_simplify(sys)
print("\nCheckpoint 6: Finished Simplifying")
prob = ODEProblem(ssys, u0, tspan, params)
print("\nCheckpoint 7: Finished creating the ODE Problem")
sol = solve(prob, Rodas4())
print("\nCheckpoint 8: Finished solving with Rodas 4")


### Timing ###
print("\nTime to convert:")
@time convert(ODESystem, complete(system))
print("\nTime to simplify:")
@time structural_simplify(sys)
print("\nTime to create the simplified problem:")
@time ODEProblem(ssys, u0, tspan, params)
print("\nTime to solve the simplified 1000 reaction system with Rodas4(): ")
@time solve(prob, Rodas4());


### Plotting ###
# C and C+
plot(sol, idxs = (0,10), lw = 3, lc = "blue", title = "Glover with UMISTrates")
plot!(sol, idxs = (0,9), lw = 3, lc = "orange", title = "Glover with UMIST rates")

#=
# CO
plot(sol, idxs = (0,16), lw = 3, lc = "green", title = "Glover")

# O 
plot(sol, idxs = (0,12), lw = 3, lc = "blue", title = "Glover")

# He+
plot(sol, idxs = (0,8), lw = 3, lc = "light pink", title = "Glover")

# CH and CH2
plot(sol, idxs = (0,17), lw = 3, lc = "blue", title = "Glover")
plot!(sol, idxs = (0,18), lw = 3, lc = "light blue", title = "Glover")

# OH, OH+, H2O, H2O+, and O2
plot(sol, idxs = (0,13), lw = 3, lc = "green", title = "Glover")
plot!(sol, idxs = (0,27), lw = 3, lc = "dark green", title = "Glover")
plot!(sol, idxs = (0,20), lw = 3, lc = "blue", title = "Glover")
plot!(sol, idxs = (0,28), lw = 3, lc = "light blue", title = "Glover")
plot!(sol, idxs = (0,21), lw = 3, lc = "orange", title = "Glover")

# H3+
plot(sol, idxs = (0,22), lw = 3, lc = "orange", title = "Glover")

# HCO+ 
plot(sol, idxs = (0,15), lw = 3, lc = "orange", title = "Glover")

# M 
plot(sol, idxs = (0,33), lw = 3, lc = "light blue", title = "Glover")
=#



# LEGEND for species ID's
# 1: H
# 2: e
# 3: H-
# 4: H2
# 5: H+
# 6: H2+
# 7: He
# 8: He+
# 9: C+
# 10: C
# 11: O+
# 12: O
# 13: OH
# 14: HOC+
# 15: HCO+
# 16: CO
# 17: CH
# 18: CH2
# 19: C2
# 20: H2O
# 21: O2
# 22: H3+
# 23: CH+
# 24: CH2+
# 25: CO+
# 26: CH3+
# 27: OH+
# 28: H2O+
# 29: H3O+
# 31: C-
# 32: O-
# 33: M