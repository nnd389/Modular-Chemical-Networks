#using Revise
#using DrWatson
using Catalyst
using DifferentialEquations
using Plots
using OrdinaryDiffEq
using ModelingToolkit
using LSODA

# TO DOs:
# FIX K9, K10, K30 - these have the n_h numbers?
# FIX K12, K17 - these have strange conditionals?
# CHECK K159 and K160 - these have k210 in their rate definitions
# FIX K165 - rate is correct, but has an H(s) in it?
# FIX K218 - has x_H2 and x_CO in the rate def, currently commented out
# check, k218 and k198 are the same reaction

# FIX av value, see equation 2 in paper
# FIX Td

#look for corrected version? look at astroph 


#@time begin
    # %% Set The timespan, parameters, and initial conditions
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs
# 30,000 yrs approximately equals 8e11 
#tspan = (0, 1000 * seconds_per_year)
#tspan = (0, 5e9)


# Nelson has CHx (CH and CH2) and OHx (OH, H2O, O2) start at zero
u0 = [  0,        # 1: H FIX
        2e-4,     # 2: e Nelson has 2e-4
        0,        # 3: H- FIX
        0.5,      # 4: H2 Nelson has 0.5
        0,        # 5: H+ FIX
        0,        # 6: H2+ FIX
        0.1,      # 7: He Nelson has 0.1
        7.866e-7, # 8: He+ Nelson has 7.866e-7
        2e-4,     # 9: C+ Nelson has 2e-4, Glover has...
        0,        # 10: C Nelson has 0, Glover has 1.41e-4, see section 4
        0,        # 11: O+ FIX
        4e-4,     # 12: O Nelson has 4e-4, Glover has 3.16e-4, see section 4
        0,        # 13: OH Nelson has 0 for OHx
        0,        # 14: HOC+ FIX
        0,        # 15: HCO+ Nelson has 0
        0,        # 16: CO Nelson has 0
        0,        # 17: CH FIX
        0,        # 18: CH2 Nelson has 0 for CHx
        0,        # 19: C2 FIX
        0,        # 20: H2O Nelson has 0
        0,        # 21: O2 Nelson has 0, OHx FIX-------CHECK IF THIS IS O2+ FIND O2
        9.059e-9, # 22: H3+ Nelson has 9.059e-9
        0,        # 23: CH+ FIX
        0,        # 24: CH2+ FIX
        0,        # 25: CO+ FIX
        0,        # 26: CH3+ Nelson has 0 for CHx
        0,        # 27: OH+ FIX
        0,        # 28: H2O+ FIX
        0,        # 29: H3O+ FIX
        0,        # 30: O2+ FIX
        0,        # 31: C- FIX
        0,        # 32: O- FIX
        2e-7      # 33: M Nelson has 2e-7
        ]



#Te = 0.06032129704818 # 700k = 0.06032129704818 ev
#Te = 0.02585198444922 # 300k = 0.02585198444922 ev
#Te = 0.008617328149741 # 100k = 0.008617328149741 ev
#Te = 0.00430866407487 # 50k = 0.00430866407487 ev
#Te = 0.0008617328149741 # 10k = 0.0008617328149741 ev






T = 10 # from glover paper, T=60 section 4, nelson has 10k
Td = 10
Te = 0.0008617328149741 # 10k = 0.0008617328149741 ev
n_H = 611 # Nelson has 611
Av = 2 # Glover Section 2 gives formula: n_H/(1.87e21)
cr_ion_rate = 6e-18 # I tried 6e-18, DESPOTIC suggests 2.0e-16









#=
# These are the parameters that make Glover looks very similar to Nelson over 30000 years
T = 110 # from glover paper, T=60 section 4, nelson has 10k
Td = 110
Te = 0.009479060964715 # 110k = 0.009479060964715 ev
n_H = 46 # Nelson has 611
Av = 2 # Glover Section 2 gives formula: n_H/(1.87e21)
cr_ion_rate = 6e-18 # I tried 6e-18, DESPOTIC suggests 2.0e-16
=#

params = Dict(
    :T => T,  
    :Te => Te,
    :Td => Td,
    :n_H => n_H, # n_H = 300, from glover paper, section 4
    :Av => Av,
    :cr_ion_rate => cr_ion_rate,
    :k1 => 10^(-17.845 + 0.762 * log10(T) + 0.1523 * (log10(T))^2 - 0.03274 * (log10(T))^3),
    :k2 => 1.5e-9,
    :k3 => 10^(-19.38 - 1.523 * log10(T) + 1.118 * (log10(T))^2 - 0.1269 * (log10(T))^3), 
    :k4 => 6.4e-10 , 
    :k5 => 2.4e-6 * (T)^(-1/2) * (1 + T/20000), 
    :k6 => 1e-8,
    :k7 => (-3.3232183e-7 + 3.3735382e-7 * log(T) - 1.4491368e-7 * (log(T))^2 + 3.4172805e-8 * (log(T))^3 - 4.781372e-9 * (log(T))^4 + 3.9731542e-10 * (log(T))^5 - 1.8171411e-11 * (log(T))^6 + 3.5311932e-13 * (log(T))^7) * exp(-21237.15/T), 
    :k8 => 3.73e-9 * (T)^(0.1121) * exp(-99430/T), 
    :k9 => 6.67e-12 * (T)^(1/2) * exp(-(1 + 63590/T)), 
    :k10 => (5.996e-30 * (T)^(4.1881)) / ((1 + 6.761e-6 * T)^(5.6881)) * exp(-54657.4/T), 
    :k11 => exp(-3.271396786e1 + 1.3536556e1 * log(Te) - 5.73932875 * (log(Te))^2 + 1.56315498 * (log(Te))^3 - 2.877056e-1 * (log(Te))^4 + 3.48255977e-2 * (log(Te))^5 - 2.63197617e-3 * (log(Te))^6 + 1.11954395e-4 * (log(Te))^7 - 2.03914985e-6 * (log(Te))^8), 
    :k12 => 1.269e-13 * ((315614/T)^(1.503)) * (1 + (604625/T)^(0.47))^(-1.923),
    :k13 => exp(-1.801849334e1 + 2.3608522 * log(Te) - 2.827443e-1 * (log(Te))^2 + 1.62331664e-2 * (log(Te))^3 - 3.36501203e-2 * (log(Te))^4 + 1.17832978e-2 * (log(Te))^5 - 1.6561947e-3 * (log(Te))^6 + 1.0682752e-4 * (log(Te))^7 - 2.63128581e-6 * (log(Te))^8), 
    :k14 => 2.5634e-9 * (Te)^(1.78186),
    :k15 => 6.9e-9 * T^(-0.35),
    :k16 => exp(-4.409864886e1 + 2.391596563e1 * log(Te) - 1.07532302e1 * (log(Te))^2 + 3.05803875 * (log(Te))^3 - 5.6851189e-1 * (log(Te))^4 + 6.79539123e-2 * (log(Te))^5 - 5.0090561e-3 * (log(Te))^6 + 2.06723616e-4 * (log(Te))^7 - 3.64916141e-6 * (log(Te))^8), 
    :k17 => 1e-11 * T^(-0.5) * (12.72 - 1.615 * log10(T) - 0.3162 * (log10(T))^2 + 0.0493 * (log10(T))^3),
    :k18 => 1.25e-15 * (T/300)^(0.25), 
    :k19 => 1.26e-9 * T^(-0.75) * exp(-127500/T),
    :k20 => 4.67e-12 * (T/300)^(-0.6), # Glover rate for C⁺ + e --> C 
    #:k20 => 2.36e-12 * (T/300)^(-0.29) * exp(-(-17.6/T)), # Umist rate 6158 for C⁺ + e --> C 
    #:k20 => 1.4e-10 * (T)^(-0.61), # Nelson rate for C⁺ + e --> C
    :k20_reverse => 3.1e-10 * exp(-(3.3*Av)), # Umist rate 5827 (Photoprocess)for C --> C⁺ + e, There are three Umist matches for this reaction, this is the correct one
    #:k20_reverse => 3e-10 * 2 * exp(-3 * Av) # Nelson rate for C --> C⁺ + e
    :k21 => 1.3e-10 * T^(-0.64),
    :k22 => 6.85e-8 * (0.193 + (11.26/Te))^(-1) * (11.26/Te)^(0.25) * exp(-11.26/Te), 
    :k23 => 3.59e-8 * (0.073 + (13.6/Te))^(-1) * (13.6/Te)^(0.34) * exp(-13.6/Te),
    :k24 => (4.99e-11 * T^(0.405) + 7.54e-10 * T^(-0.458)),
    :k25 => (1.08e-11 * T^(0.517) + 4e-10 * T^(0.00669)) * exp(-227/T),
    :k26 => (4.991e-15 * (T/10000)^(0.3794) * exp(-T/1121000) + 2.78e-15 * (T/10000)^(-0.2163)) * exp(T/815800),
    :k27 => 3.9e-16 * T^(0.213),
    :k28 => 6.08e-14 * (T/10000)^(1.96) * exp(-170000/T),
    :k29 => 8.58e-17 * (T)^(0.757),
    :k30 => 10^(-27.029 + 3.801 * log10(T) - 29487/T), 
    :k31 => 6e-9 * exp(-50900/T), 
    :k32 => 3.8e-10, 
    :k33 => 4e-10, 
    :k34 => 6.64e-10 * exp(-11700/T), 
    :k35 => 1.31e-10 * exp(-80/T), 
    :k36 => 5.46e-10 * exp(-1943/T), 
    :k37 => 6.59e-11, 
    :k38 => 6.6e-11,
    :k39 => 6.64e-11, 
    :k40 => 1.33e-10, 
    :k41 => 8e-11, 
    :k42 => 5e-11 * (T/300)^(0.5),
    :k43 => 3.14e-13 * (T/300)^(2.7) * exp(-3150/T), 
    :k44 => 6.99e-14 * (T/300)^(2.8) * exp(-1950/T), 
    :k45 => 2.05e-12 * (T/300)^(1.52) * exp(-1736/T), 
    :k46 => 1e-10, 
    :k47 => 3.5e-11,
    :k48 => 1.65e-12 * (T/300)^(1.14) * exp(-50/T), 
    :k49 => 1.59e-11 * (T/300)^(1.2) * exp(-9610/T), 
    :k50 => 2.61e-10 * exp(-8156/T), 
    :k51 => 3.16e-10 * exp(-21890/T), 
    :k52 => 4.7e-11 * (T/300)^(-0.34),
    :k53 => 1.1e-10 * (T/300)^(0.5) * exp(-77700/T) , 
    :k54 => 2.24e-9 * (T/300)^(0.042) * exp(-T/46600), 
    :k55 => 7.7e-9 * exp(-17560/T), 
    :k56 => 2.4e-9, 
    :k57 => 2e-9, 
    :k58 => 1e-10 * exp(-4640/T), 
    :k59 => 7.5e-10, 
    :k60 => 1.2e-9, 
    :k61 => 3.5e-10, 
    :k62 => 1.4e-9, 
    :k63 => 1e-9 * exp(-7080/T), 
    :k64 => 1.6e-9, 
    :k65 => 7.5e-10, # change back to -9?
    :k66 => 7e-10 * exp(-10560/T), 
    :k67 => 4e-10, 
    :k68 => 4.8e-10, 
    :k69 => 1.7e-9, 
    :k70 => 1.5e-9 , 
    :k71 => 8.4e-10, 
    :k72 => 1.3e-9, 
    :k73 => 7.7e-10, 
    :k74 => 1.01e-9, 
    :k75 => 6.4e-10, 
    :k76 => 5.9e-9, 
    :k77 => 9e-10 , 
    :k78 => 1.8e-9, 
    :k79 => 1e-11, 
    :k80 => 3.8e-10, 
    :k81 => 6.2e-10, 
    :k82 => 9.1e-10, 
    :k83 => 5.2e-11, 
    :k84 => 2.7e-11, 
    :k85 => 1.7e-9 , 
    :k86 => 1.1e-9, 
    :k87 => 2.5e-9, 
    :k88 => 7.2e-15, 
    :k89 => 3.7e-14 * exp(-35/T), 
    :k90 => 1.9e-9, 
    :k91 => 1.4e-9, 
    :k92 => 7.5e-10, 
    :k93 => 1.6e-9, 
    :k94 => 2.1e-9, 
    :k95 => 1.1e-9, 
    :k96 => 6.9e-9, 
    :k97 => 2.04e-10, 
    :k98 => 2.86e-10, 
    :k99 => 6.05e-11, 
    :k100 => 2e-9, 
    :k101 => 3.3e-11, 
    :k102 => 1.1e-9, 
    :k103 => 5.2e-11, 
    :k104 => 1.4e-9*(T/300)^(-0.5), 
    :k105 => 1.4e-16 * (T/300)^(-0.5), 
    :k106 => 7.5e-10, 
    :k107 => 2.3e-7 * (T/300)^(-0.5), 
    :k108 => 2.3e-7 * (T/300)^(-0.5), 
    :k109 => 2.32e-7 * (T/300)^(-0.52) * exp(T/22400), 
    :k110 => 2.34e-8 * (T/300)^(-0.52), 
    :k111 => 4.36e-8 * (T/300)^(-0.52), 
    :k112 => 7e-8 * (T/300)^(-0.5), 
    :k113 => 1.6e-7 * (T/300)^(-0.6), 
    :k114 => 4.03e-7 * (T/300)^(-0.6), 
    :k115 => 7.68e-8 * (T/300)^(-0.6), 
    :k116 => 7.75e-8 * (T/300)^(-0.5), 
    :k117 => 1.95e-7 * (T/300)^(-0.5), 
    :k118 => 2e-7 * (T/300)^(-0.4), 
    :k119 => 6.3e-9 * (T/300)^(-0.48), 
    :k120 => 3.05e-7 * (T/300)^(-0.5), 
    :k121 => 3.9e-8 * (T/300)^(-0.5), 
    :k122 => 8.6e-8 * (T/300)^(-0.5), 
    :k123 => 1.08e-7 * (T/300)^(-0.5), 
    :k124 => 6.02e-8 * (T/300)^(-0.5), 
    :k125 => 2.58e-7 * (T/300)^(-0.5), 
    :k126 => 5.6e-9 * (T/300)^(-0.5), 
    :k127 => 1.95e-7 * (T/300)^(-0.7), 
    :k128 => 2.75e-7 * (T/300)^(-0.55), 
    :k129 => 2.76e-7 * (T/300)^(-0.64), 
    :k130 => 2.4e-8 * (T/300)^(-0.64), 
    :k131 => 1.1e-7 * (T/300)^(-1), 
    :k132 => 1e-9, 
    :k133 => 1e-9, 
    :k134 => 1e-10, 
    :k135 => 5e-10, 
    :k136 => 1e-13, 
    :k137 => 5e-10, 
    :k138 => 5e-10, 
    :k139 => 7e-10 , 
    :k140 => 5e-10, 
    :k141 => 1e-16, 
    :k142 => 2.25e-15, 
    :k143 => 1e-17, 
    :k144 => 1e-17, 
    :k145 => 4.36e-18 * (T/300)^(0.35) * exp(-161.3/T), 
    :k146 => 2.1e-19,
    :k147 => 4.46e-16 * T^(-0.5) * exp(-4.93/(T^(2/3))), 
    :k148 => 4e-16 *(T/300)^(-0.2), 
    :k149 => 2.5e-18,
    :k150 => 1.5e-15, 
    :k151 => 9.9e-19 * (T/300)^(-0.38), 
    :k152 => 4.9e-20 * (T/300)^1.58, 
    :k153 => 5.26e-18 *(T/300)^(-5.22) * exp(-90/T), 
    :k154 => 1.32e-32 * (T/300)^-0.38,
    :k155 => 2.8e-31 * (T)^-0.6,
    :k156 => 6.9e-32 * (T)^(-0.4) ,
    :k157 => 5.99e-33 * (T/5000)^(-1.6), 
    :k158 => 6.16e-29 * (T/300)^(-3.08),
    :k159 => 100 * cr_ion_rate * 960, 
    :k160 => 100 * cr_ion_rate * 960, 
    :k161 => 4.33e-32 * (T/300)^(-1), 
    :k162 => 2.56e-31 * (T/300)^(-2), 
    :k163 => 9.2e-34 * (T/300)^(-1), 
    :k164 => 2e-11 * (T/300)^(0.44), 
    :k165 => 3e-18 * (T)^(0.5) * (1 + 1e4 * exp(-600/Td))^(-1) * (1 + 0.04*(T + Td)^(0.5) + 0.002*T + 8e-6 * T^2)^(-1), 
    :R166 => 7.1e-7 * exp(-0.5 * Av), 
    :R167 => 1.1e-9 * exp(-1.9 * Av), 
    :R168 => 5.6e-11 * exp(-2.55*Av + 0.0165*(Av)^2), 
    :R169 => 4.9e-13 * exp(-1.8 * Av), 
    :R170 => 4.9e-13 * exp(-2.3 * Av), 
    :R171 => 3.1e-10 * exp(-3.0 * Av), 
    :R172 => 2.4e-7 * exp(-0.9 * Av), 
    :R173 => 8.7e-10 * exp(-1.2 * Av), 
    :R174 => 7.7e-10 * exp(-2.8 * Av), 
    :R175 => 2.6e-10 * exp(-2.5 * Av), 
    :R176 => 7.1e-10 * exp(-1.7 * Av), 
    :R177 => 5.9e-10 * exp(-2.3 * Av), 
    :R178 => 4.6e-10 * exp(-1.7 * Av), 
    :R179 => 1.0e-9 * exp(-1.7 * Av), 
    :R180 => 1.0e-9 * exp(-1.7 * Av), 
    :R181 => 1.5e-10 * exp(-2.1 * Av), 
    :R182 => 2.4e-7 * exp(-0.5 * Av), 
    :R183 => 3.7e-10 * exp(-1.7 * Av), 
    :R184 => 1.6e-12 * exp(-3.1 * Av), 
    :R185 => 1.0e-12 * exp(-1.8 * Av), 
    :R186 => 6.0e-10 * exp(-1.7 * Av), 
    :R187 => 3.2e-11 * exp(-3.9 * Av), 
    :R188 => 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R189 => 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R190 => 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R191 => 1.5e-10 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R192 => 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R193 => 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R194 => 7.5e-12 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R195 => 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2), 
    :R196 => 5.6e-11 * exp(-3.7 * Av), 
    :R197 => 7.0e-10 * exp(-1.8 * Av), 
    :R198 => 2e-10 * exp(-3.5 * Av), 
    :R199 => cr_ion_rate * 1, 
    :R200 => cr_ion_rate * 1.1333333, 
    :R201 => cr_ion_rate * 0.037, 
    :R202 => cr_ion_rate * 0.22, 
    :R203 => cr_ion_rate * 6.5e-4, 
    :R204 => cr_ion_rate * 2, 
    :R205 => cr_ion_rate * 3.8, 
    :R206 => cr_ion_rate * 5.7, 
    :R207 => cr_ion_rate * 6.5, 
    :R208 => cr_ion_rate * 2800, 
    :R209 => cr_ion_rate * 4000, 
    :R210 => cr_ion_rate * 960, 
    :R211 => cr_ion_rate * 2700, 
    :R212 => cr_ion_rate * 2700, 
    :R213 => cr_ion_rate * 1300, 
    :R214 => cr_ion_rate * 2800, 
    :R215 => cr_ion_rate * 5300, 
    :R216 => cr_ion_rate * 4100, 
    :R217 => cr_ion_rate * 640 
    #:R218 => , 
    )


if T > 6000
    params[:k1] = 10^(-16.420 + 0.1998 * (log10(T))^2 - 5.447e-3 * (log10(T))^4 + 4.0415e-5 * (log10(T))^6)
end

if T > 300
    params[:k2] = 4e-9 * T^(-0.17)
end

if T > 617
    params[:k6] = 1.32e-6 * (T)^(-0.76)
end

if Te <= 0.1
    params[:k14] = exp(-2.0372609e1 + 1.13944933 * log(Te) - 1.4210135e-1 * (log(Te))^2 + 8.4644554e-3 * (log(Te))^3 - 1.4327641e-3 * (log(Te))^4 + 2.0122503e-4 * (log(Te))^5 + 8.6639632e-5 * (log(Te))^6 - 2.5850097e-5 * (log(Te))^7 +2.4555012e-6 * (log(Te))^8 - 8.0683825e-8 * (log(Te))^9)
end

if T > 8000
    params[:k15] = 9.6e-7*T^(-0.9)
end

if T > 10000
    params[:k19] = 4e-37 * (T)^(4.74)
end

if T > 7950 && T <= 21140
    params[:k20] = 1.23e-17 * (T/300)^(2.49) * exp(21845.6/T)
end

if T > 21140
    params[:k20] = 9.62e-8 * (T/300)^(-1.37) * exp(-115786.2/T)
end

if T > 400
    params[:k21] = 1.41e-10 * (T)^(-0.66) + 7.4e-4 * T^(-1.5)*exp(-175000/T) * (1+0.062*exp(-145000/T))
end

if T > 200 && T <= 2000
    params[:k29] = 3.25e-17 * T^(0.968)
end

if T > 2000
    params[:k29] = 2.77e-19 * T^(1.597)
end

if T > 2000
    params[:k38] = 1.02e-10 * exp(-914/T)
end

if T > 300
    params[:k42] = 5e-11 * (T/300)^(0.757) 
end

if T > 261
    params[:k47] = 1.77e-11 * exp(178/T)
end

if T > 295
    params[:k52] = 2.48e-12 * (T/300)^(1.54)*exp(613/T)
end

if T > 300
    params[:k146] = 3.09e-17 * (T/300)^(0.33) * exp(-1629/T)
end

if T > 300
    params[:k149] = 3.14e-18 *(T/300)^(-0.15) * exp(68/T)
end

if T > 300
    params[:k154] = 1.32e-32 * (T/300)^-1
end

if T > 5000
    params[:k157] = 5.99e-33 * (T/5000)^(-0.64) * exp(5255/T)
end

if T > 2000
    params[:k158] = 2.14e-29 * (T/300)^(-3.08) * exp(2114/T)
end

if Av > 15
    params[:R188] = 5.0e-11 * exp(-2.8 * Av)
    params[:R189] = 5.0e-11 * exp(-2.8 * Av)
    params[:R190] = 5.0e-10 * exp(-2.8 * Av)
    params[:R191] = 1.5e-10 * exp(-2.8 * Av)
    params[:R192] = 2.5e-11 * exp(-2.8 * Av)
    params[:R193] = 2.5e-11 * exp(-2.8 * Av)
    params[:R194] = 7.5e-12 * exp(-2.8 * Av)
    params[:R195] = 2.5e-11 * exp(-2.8 * Av)
end



# %% Network admin things
#allvars = @strdict tspan u0 params
#allvars = @strdict tspan u0

@variables t 
@species H⁻(t) H2⁺(t) H3⁺(t) CH⁺(t) CH2⁺(t) OH⁺(t) H2O⁺(t) H3O⁺(t) CO⁺(t) HOC⁺(t) O⁻(t) C⁻(t) O2⁺(t) e(t) H⁺(t) H(t) H2(t) He(t) He⁺(t) C(t) C⁺(t) O(t) O⁺(t) OH(t) H2O(t) CO(t) C2(t) O2(t) H2O⁺(t) CH(t) CH2(t) CH3⁺(t) M(t)
@parameters T n_H Av k1 k2 k6 k12 k14 k15 k17 k19 k20 k21 k29 k38 k42 k47 k52 k146 k149 k154 k157 k158 R188 R189 R190 R191 R192 R193 R194 R195


# %% Define the network
reaction_equations = [
    (@reaction n_H * k1, H + e --> H⁻),    #k1
    (@reaction n_H * k2, H⁻ + H --> H2 + e), #k2
    (@reaction n_H * k3, H + H⁺ --> H2⁺ ), 
    (@reaction n_H * k4, H + H2⁺ --> H2 + H⁺), 
    (@reaction n_H * k5, H⁻ + H⁺--> H + H), 
    (@reaction n_H * k6, H2⁺ + e --> H + H), 
    (@reaction n_H * k7, H2 + H⁺ --> H2⁺ + H), 
    (@reaction n_H * k8, H2 + e --> H + H + e), 
    (@reaction n_H * k9, H2 + H --> H + H + H), 
    (@reaction n_H * k10, H2 + H2 --> H2 + H + H), 
    (@reaction n_H * k11, H + e --> H⁺ + e + e), 
    (@reaction n_H * k12, H⁺ + e --> H ),
    (@reaction n_H * k13, H⁻ + e --> H + e + e), 
    (@reaction n_H * k14, H⁻ + H --> H + H + e), 
    (@reaction n_H * k15, H⁻ + H⁺--> H2⁺ + e),
    (@reaction n_H * k16, He + e --> He⁺ + e + e), 
    (@reaction n_H * k17, He⁺ + e --> He ), 
    (@reaction n_H * k18, He⁺ + H --> He + H⁺), 
    (@reaction n_H * k19, He + H⁺--> He⁺ + H), 
    (@reaction n_H * k20, C⁺ + e --> C ), # See parameters to switch out Glover rate for Nelson or Umist rate
    (@reaction k20_reverse, C --> C⁺ + e ), # Umist rate 5827 (Photoprocess)for C --> C⁺ + e, There are three Umist matches for this reaction, this is the correct one. See parameters to switch out for Nelson rate. 
    (@reaction n_H * k21, O⁺ + e --> O ), 
    (@reaction n_H * k22, C + e --> C⁺ + e + e), 
    (@reaction n_H * k23, O + e --> O⁺ + e + e), 
    (@reaction n_H * k24, O⁺ + H --> O + H⁺), 
    (@reaction n_H * k25, O + H⁺--> O⁺ + H), 
    (@reaction n_H * k26, O + He⁺ --> O⁺ + He), 
    (@reaction n_H * k27, C + H⁺ --> C⁺ + H), 
    (@reaction n_H * k28, C⁺ + H --> C + H⁺), 
    (@reaction n_H * k29, C + He⁺ --> C⁺ + He), 
    (@reaction n_H * k30, H2 + He --> H + H + He), 
    (@reaction n_H * k31, OH + H --> O + H + H), 
    (@reaction n_H * k32, HOC⁺ + H2 --> HCO⁺ + H2), 
    (@reaction n_H * k33, HOC⁺ + CO --> HCO⁺ + CO), 
    (@reaction n_H * k34, C + H2 --> CH + H), 
    (@reaction n_H * k35, CH + H --> C + H2), 
    (@reaction n_H * k36, CH + H2 --> CH2 + H), 
    (@reaction n_H * k37, CH + C --> C2 + H), 
    (@reaction n_H * k38, CH+O --> CO+H), # Nelson Match Neutral-Neutral reaction 1, Nelson has 2e-10
    (@reaction n_H * k39, CH2 + H --> CH + H2), 
    (@reaction n_H * k40, CH2 + O --> CO + H + H), 
    (@reaction n_H * k41, CH2 + O --> CO + H2), 
    (@reaction n_H * k42, C2 + O --> CO + C), 
    (@reaction n_H * k43, O + H2 --> OH + H), 
    (@reaction n_H * k44, OH + H --> O + H2), 
    (@reaction n_H * k45, OH + H2 --> H2O + H), 
    (@reaction n_H * k46, OH + C --> CO + H), # Nelson Match Neutral-Neutral reaction 2, Nelson has 5.8e-12 * T^(0.5)
    (@reaction n_H * k47, OH + O --> O2 + H), 
    (@reaction n_H * k48, OH + OH --> H2O + H), 
    (@reaction n_H * k49, H2O + H --> H2 + OH), 
    (@reaction n_H * k50, O2 + H --> OH + O), 
    (@reaction n_H * k51, O2 + H2 --> OH + OH), 
    (@reaction n_H * k52, O2 + C --> CO + O), 
    (@reaction n_H * k53, CO + H --> C + OH), 
    (@reaction n_H * k54, H2⁺ + H2 --> H3⁺ + H), 
    (@reaction n_H * k55, H3⁺ + H --> H2⁺ + H2), 
    (@reaction n_H * k56, C + H2⁺ --> CH⁺ + H), 
    (@reaction n_H * k57, C + H3⁺ --> CH⁺ + H2), # Nelson Match Ion-Molecule reaction 1
    (@reaction n_H * k58, C⁺ + H2 --> CH⁺ + H), # Nelson Match Ion-Molecule reaction 6, Nelson has 4e-16
    (@reaction n_H * k59, CH⁺ + H --> C⁺ + H2), 
    (@reaction n_H * k60, CH⁺ + H2 --> CH2⁺ + H), 
    (@reaction n_H * k61, CH⁺ + O --> CO⁺ + H), 
    (@reaction n_H * k62, CH2 + H⁺--> CH⁺ + H2), 
    (@reaction n_H * k63, CH2⁺ + H --> CH⁺ + H2), 
    (@reaction n_H * k64, CH2⁺ + H2 --> CH3⁺ + H), 
    (@reaction n_H * k65, CH2⁺ + O --> HCO⁺ + H), 
    (@reaction n_H * k66, CH3⁺ + H --> CH2⁺ + H2), 
    (@reaction n_H * k67, CH3⁺ + O --> HCO⁺ + H2), 
    (@reaction n_H * k68, C2 + O⁺ --> CO⁺ + C), 
    (@reaction n_H * k69, O⁺ + H2 --> OH⁺ + H), 
    (@reaction n_H * k70, O + H2⁺ --> OH⁺ + H), 
    (@reaction n_H * k71, O + H3⁺ --> OH⁺ + H2), # Nelson Match Ion-Molecule reaction 2 Nelson has 8e-10
    (@reaction n_H * k72, OH + H3⁺ --> H2O⁺ + H2), 
    (@reaction n_H * k73, OH + C⁺ --> CO⁺ + H), 
    (@reaction n_H * k74, OH⁺ + H2 --> H2O⁺ + H), 
    (@reaction n_H * k75, H2O⁺ + H2 --> H3O⁺ + H), 
    (@reaction n_H * k76, H2O + H3⁺ --> H3O⁺ + H2), 
    (@reaction n_H * k77, H2O + C⁺ --> HCO⁺ + H), 
    (@reaction n_H * k78, H2O + C⁺ --> HOC⁺ + H), 
    (@reaction n_H * k79, H3O⁺ + C --> HCO⁺ + H2), 
    (@reaction n_H * k80, O2 + C⁺ --> CO⁺ + O), 
    (@reaction n_H * k81, O2 + C⁺ --> CO + O⁺), 
    (@reaction n_H * k82, O2 + CH2⁺ --> HCO⁺ + OH), 
    (@reaction n_H * k83, O2⁺ + C --> CO⁺ + O), 
    (@reaction n_H * k84, CO + H3⁺ --> HOC⁺ + H2), 
    (@reaction n_H * k85, CO + H3⁺ --> HCO⁺ + H2), # Nelson Match Ion-Molecule reaction 3
    (@reaction n_H * k86, HCO⁺ + C --> CO + CH⁺), 
    (@reaction n_H * k87, HCO⁺ + H2O --> CO + H3O⁺), 
    (@reaction n_H * k88, H2 + He⁺  -->  He + H2⁺), # Nelson Match Ion-Molecule reaction 4, Nelson has 7e-15
    (@reaction n_H * k89, H2 + He⁺ --> He + H + H⁺), # Nelson Match Ion-Molecule reaction 4, Nelson has 7e-15
    (@reaction n_H * k90, CH + H⁺ -->  CH⁺ + H), 
    (@reaction n_H * k91, CH2 + H⁺ -->  CH2⁺ + H), 
    (@reaction n_H * k92, CH2 + He⁺ -->  C⁺ + He + H2), 
    (@reaction n_H * k93, C2 + He⁺ -->  C⁺ + C + He), 
    (@reaction n_H * k94, OH + H⁺ -->  OH⁺ + H), 
    (@reaction n_H * k95, OH + He⁺ -->  O⁺ + He + H), 
    (@reaction n_H * k96, H2O + H⁺ -->  H2O⁺ + H), 
    (@reaction n_H * k97, H2O + He⁺ -->  OH + He + H⁺), 
    (@reaction n_H * k98, H2O + He⁺ -->  OH⁺ + He + H), 
    (@reaction n_H * k99, H2O + He⁺ -->  H2O⁺ + He), 
    (@reaction n_H * k100, O2 + H⁺ -->  O2⁺ + H), 
    (@reaction n_H * k101, O2 + He⁺ -->  O2⁺ + He), 
    (@reaction n_H * k102, O2 + He⁺ -->  O⁺ + O + He), 
    (@reaction n_H * k103, O2⁺ + C  -->  O2 + C⁺), 
    (@reaction n_H * k104, CO + He⁺ -->  C⁺ + O + He), # Nelson Match Ion-Molecule reaction 5, Nelson has 1.6e-9
    (@reaction n_H * k105, CO + He⁺ -->  C + O⁺ + He), 
    (@reaction n_H * k106, CO⁺ + H  -->  CO + H⁺), 
    (@reaction n_H * k107, C⁻ + H⁺ -->  C + H), 
    (@reaction n_H * k108, O⁻ + H⁺ -->  O + H), 
    (@reaction n_H * k109, He⁺ + H⁻  -->  He + H), 
    (@reaction n_H * k110, H3⁺ + e --> H2 + H), 
    (@reaction n_H * k111, H3⁺ + e --> H + H + H), 
    (@reaction n_H * k112, CH⁺ + e --> C + H), 
    (@reaction n_H * k113, CH2⁺ + e --> CH + H), 
    (@reaction n_H * k114, CH2⁺ + e --> C + H + H), 
    (@reaction n_H * k115, CH2⁺ + e --> C + H2), 
    (@reaction n_H * k116, CH3⁺ + e --> CH2 + H), 
    (@reaction n_H * k117, CH3⁺ + e --> CH + H2), 
    (@reaction n_H * k118, CH3⁺ + e --> CH + H + H), 
    (@reaction n_H * k119, OH⁺ + e --> O + H), 
    (@reaction n_H * k120, H2O⁺ + e --> O + H + H), 
    (@reaction n_H * k121, H2O⁺ + e --> O+H2), 
    (@reaction n_H * k122, H2O⁺ + e --> OH + H), 
    (@reaction n_H * k123, H3O⁺ + e --> H + H2O), 
    (@reaction n_H * k124, H3O⁺ + e --> OH + H2), 
    (@reaction n_H * k125, H3O⁺ + e --> OH + H + H), 
    (@reaction n_H * k126, H3O⁺ + e --> O + H + H2), 
    (@reaction n_H * k127, O2⁺ + e --> O + O), 
    (@reaction n_H * k128, CO⁺ + e --> C + O), 
    (@reaction n_H * k129, HCO⁺ + e --> CO + H), 
    (@reaction n_H * k130, HCO⁺ + e --> OH + C), 
    (@reaction n_H * k131, HOC⁺ + e --> CO + H), 
    (@reaction n_H * k132, H⁻ + C --> CH + e), 
    (@reaction n_H * k133, H⁻ + O --> OH + e), 
    (@reaction n_H * k134, H⁻ + OH --> H2O + e), 
    (@reaction n_H * k135, C⁻ + H --> CH + e), 
    (@reaction n_H * k136, C⁻ + H2 --> CH2 + e), 
    (@reaction n_H * k137, C⁻ + O -->  CO + e), 
    (@reaction n_H * k138, O⁻ + H  -->  OH + e), 
    (@reaction n_H * k139, O⁻ + H2 -->  H2O + e), 
    (@reaction n_H * k140, O⁻ + C --> CO + e), 
    (@reaction n_H * k141, H2 + H⁺ --> H3⁺), 
    (@reaction n_H * k142, C + e --> C⁻ ), 
    (@reaction n_H * k143, C + H --> CH ), 
    (@reaction n_H * k144, C + H2 --> CH2 ), 
    (@reaction n_H * k145, C + C --> C2 ), 
    (@reaction n_H * k146, C + O  --> CO ), 
    (@reaction n_H * k147, C⁺ + H --> CH⁺ ), 
    (@reaction n_H * k148, C⁺ + H2 --> CH2⁺ ), 
    (@reaction n_H * k149, C⁺ + O --> CO⁺ ), 
    (@reaction n_H * k150, O + e --> O⁻ ), 
    (@reaction n_H * k151, O+H --> OH ), 
    (@reaction n_H * k152, O + O --> O2 ), 
    (@reaction n_H * k153, OH + H --> H2O ), 
    (@reaction n_H^2 * k154, H + H + H --> H2 + H), 
    (@reaction n_H^2 * k155, H + H + H2 --> H2 + H2), 
    (@reaction n_H^2 * k156, H + H + He --> H2 + He), 
    (@reaction n_H^2 * k157, C + C + M --> C2 + M), 
    (@reaction n_H^2 * k158, C + O + M --> CO + M), 
    (@reaction n_H^2 * k159, C⁺ + O + M --> CO⁺ + M), # FIX: these rely on k210, which doesn't exist?
    (@reaction n_H^2 * k160, C + O⁺ + M --> CO⁺ + M), # FIX: these rely on k210, which doesn't exist?
    (@reaction n_H^2 * k161, O + H + M --> OH + M), 
    (@reaction n_H^2 * k162, OH + H + M --> H2O + M), 
    (@reaction n_H^2 * k163, O + O + M --> O2 + M), 
    (@reaction n_H * k164, O + CH -->  HCO⁺ + e), 
    (@reaction n_H * k165, H+H --> H2), # FIX H is H(s) here?
    


    # Photochemical Reactions (R166) start below
    (@reaction R166, H⁻ --> H + e), 
    (@reaction R167, H2⁺ --> H + H⁺), 
    (@reaction R168, H2   -->  H+H), 
    (@reaction R169, H3⁺ --> H2 + H⁺),
    (@reaction R170, H3⁺ --> H2⁺ + H), 
    (@reaction R171, C -->  C⁺ + e), #171
    (@reaction R172, C⁻ --> C + e), 
    (@reaction R173, CH --> C + H), 
    (@reaction R174, CH --> CH⁺ + e), 
    (@reaction R175, CH⁺ --> C + H⁺), 
    (@reaction R176, CH2 --> CH + H), 
    (@reaction R177, CH2 --> CH2⁺ + e), 
    (@reaction R178, CH2⁺ --> CH⁺ + H), 
    (@reaction R179, CH3⁺ --> CH2⁺ + H), 
    (@reaction R180, CH3⁺ --> CH⁺ + H2), 
    (@reaction R181, C2 --> C + C), 
    (@reaction R182, O⁻ --> O + e),
    (@reaction R183, OH --> O + H), 
    (@reaction R184, OH --> OH⁺ + e), 
    (@reaction R185, OH⁺  --> O + H⁺), 
    (@reaction R186, H2O --> OH + H), 
    (@reaction R187, H2O --> H2O⁺ + e), 
    (@reaction R188, H2O⁺ --> H2⁺ + O), 
    (@reaction R189, H2O⁺ --> H⁺ + OH), 
    (@reaction R190, H2O⁺ --> O⁺ + H2), 
    (@reaction R191, H2O⁺ --> OH⁺ + H), 
    (@reaction R192, H3O⁺ --> H⁺ + H2O), 
    (@reaction R193, H3O⁺ --> H2⁺ + OH), 
    (@reaction R194, H3O⁺ --> H2O⁺ + H), 
    (@reaction R195, H3O⁺ --> OH⁺ + H2), 
    (@reaction R196, O2 --> O2⁺ + e), 
    (@reaction R197, O2 --> O + O), 
    (@reaction R198, CO   --> C + O), #from UMIST


    # Cosmic Ray Reactions (R199) start below 
    (@reaction R199, H --> H⁺ + e), 
    (@reaction R200, He --> He⁺ + e), # cr_rate = 6e-18, Nelson match CR Ionization reaction 2
    (@reaction R201, H2 --> H⁺ + H + e), 
    (@reaction R202, H2 --> H + H), 
    (@reaction R203, H2 --> H⁺ + H⁻), 
    (@reaction R204, H2 --> H2⁺ + e),  # cr_rate = 6e-18, Nelson match CR Ionization reaction 1
    (@reaction R205, C --> C⁺ + e), 
    (@reaction R206, O --> O⁺ + e), 
    (@reaction R207, CO --> CO⁺ + e), 
    (@reaction R208, C --> C⁺ + e), 
    (@reaction R209, CH --> C + H), 
    (@reaction R210, CH⁺ --> C⁺ + H), 
    (@reaction R211, CH2 --> CH2⁺ + e), 
    (@reaction R212, CH2 --> CH + H), 
    (@reaction R213, C2 --> C + C), 
    (@reaction R214, OH --> O + H), 
    (@reaction R215, H2O --> OH + H), 
    (@reaction R216, O2 --> O + O), 
    (@reaction R217, O2 --> O2⁺ + e), 
    #(@reaction cr_ion_rate * 2e-10 * exp(-3.5 * Av), CO   --> C + O) # this is a cosmic ray reaction, NOT photoreaction
    #(@reaction 10e-10 * (1) * (1.7) * exp(-3 *  Av) , CO --> C + O) # I copied the rate for this one from nelson
]


### Turn the Network into a system of ODEs ###
@named system = ReactionSystem(reaction_equations, t)
print("\nCheckpoint 4: Finished creating the reaction system")
sys = convert(ODESystem, complete(system))
print("\nCheckpoint 5: Finished converting to an ODE System")

simplified_sys = structural_simplify(sys)
completed_sys = complete(sys)
print("\nCheckpoint 6: Finished Simplifying")

prob_simplify = ODEProblem(simplified_sys, u0, tspan, params)
prob_complete = ODEProblem(completed_sys, u0, tspan, params)
print("\nCheckpoint 7: Finished creating the ODE Problem")

#sol = solve(prob, Rodas4(), saveat=1e11)
sol_simplify = solve(prob_simplify, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e9)
sol_complete = solve(prob_complete, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e9)

#sol = solve(prob, lsoda(), saveat=1e10)
#print("\nCheckpoint 8: Finished solving with Rodas 4")
print("\n")

### Timing ###
print("\nGlover:")
print("\nTime to convert:")
@time convert(ODESystem, complete(system))

print("\nTime to simplify:")
@time structural_simplify(sys)
print("Time to complete:")
@time complete(sys)

print("\nTime to create the simplified problem:")
@time ODEProblem(simplified_sys, u0, tspan, params)
print("Time to create the completed problem:")
@time ODEProblem(completed_sys, u0, tspan, params)

print("\nTime to solve the simplified system with Rodas4(): ")
@time solve(prob_simplify, Rodas4());
print("Time to solve the completed system with Rodas4(): ")
@time solve(prob_complete, Rodas4());


plot(sol_complete, idxs = (0,10), lw = 3, lc = "blue")
plot!(sol_complete, idxs = (0,9), lw = 3, lc = "orange", title = "Glover with Glover rates: C and C+")


#=
### Plotting ###
# C and C+
plot(sol, idxs = (0,10), lw = 3, lc = "blue")
plot!(sol, idxs = (0,9), lw = 3, lc = "orange", title = "Glover with Glover rates: C and C+")
savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/C_and_Cp_Glover.png")

# CO
plot(sol, idxs = (0,16), lw = 3, lc = "green", title = "Glover with Glover rates: CO")
savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/CO_Glover.png")

# O 
plot(sol, idxs = (0,12), lw = 3, lc = "blue", title = "Glover with Glover rates: O")
savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/O_Glover.png")

# He+
plot(sol, idxs = (0,8), lw = 3, lc = "light pink", title = "Glover with Glover rates: He+")
savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/Hep_Glover.png")

# CH and CH2 = CHx
plot(sol, idxs = (0,17), lw = 3, lc = "blue")
plot!(sol, idxs = (0,18), lw = 3, lc = "light blue", title = "Glover with Glover rates: CHx")
savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/CH_and_CH2_Glover.png")

# OH, OH+, H2O, H2O+, and O2 = OHx
plot(sol, idxs = (0,13), lw = 3, lc = "green")
plot!(sol, idxs = (0,27), lw = 3, lc = "dark green")
plot!(sol, idxs = (0,20), lw = 3, lc = "blue")
plot!(sol, idxs = (0,28), lw = 3, lc = "light blue")
plot!(sol, idxs = (0,21), lw = 3, lc = "orange", title = "Glover with Glover rates: OHx")
savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/OH_OHp_H2O_H2Op_O2_Glover.png")

# H3+
plot(sol, idxs = (0,22), lw = 3, lc = "orange", title = "Glover with Glover rates: H3+")
savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/H3p_Glover.png")

# HCO+ 
plot(sol, idxs = (0,15), lw = 3, lc = "orange", title = "Glover with Glover rates: HCO+")
savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/HCOp_Glover.png")

# M 
plot(sol, idxs = (0,33), lw = 3, lc = "light blue", title = "Glover with Glover rates: M")
savefig("/Users/kneenaugh/Desktop/Git/AstroChemNetwork/plots/Mg_Fe_Na_Glover.png")

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





