using Revise
using DrWatson
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

print("Checkpoint 1\n")
#@time begin
    # %% Set The timespan, parameters, and initial conditions
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs
# 30,000 yrs approximately equals 8e11 
#tspan = (0, 1000 * seconds_per_year)
#tspan = (0, 5e9)
print("Checkpoint 2\n")

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
    :k6 => 1e-8,
    :k12 => 1.269e-13 * ((315614/T)^(1.503)) * (1 + (604625/T)^(0.47))^(-1.923),
    :k14 => 2.5634e-9 * (Te)^(1.78186),
    :k15 => 6.9e-9 * T^(-0.35),
    :k17 => 1e-11 * T^(-0.5) * (12.72 - 1.615 * log10(T) - 0.3162 * (log10(T))^2 + 0.0493 * (log10(T))^3),
    :k19 => 1.26e-9 * T^(-0.75) * exp(-127500/T),
    :k20 => 4.67e-12 * (T/300)^(-0.6),
    :k21 => 1.3e-10 * T^(-0.64),
    :k29 => 8.58e-17 * (T)^(0.757),
    :k38 => 6.6e-11,
    :k42 => 5e-11 * (T/300)^(0.5),
    :k47 => 3.5e-11,
    :k52 => 4.7e-11 * (T/300)^(-0.34),
    :k146 => 2.1e-19,
    :k149 => 2.5e-18,
    :k154 => 1.32e-32 * (T/300)^-0.38,
    :k157 => 5.99e-33 * (T/5000)^(-1.6), 
    :k158 => 6.16e-29 * (T/300)^(-3.08),
    :R188 => 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R189 => 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R190 => 5.0e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R191 => 1.5e-10 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R192 => 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R193 => 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R194 => 7.5e-12 * exp(-2.55*Av + 0.0165*(Av)^2),
    :R195 => 2.5e-11 * exp(-2.55*Av + 0.0165*(Av)^2)
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
print("Checkpoint 3\n")
@variables t 
@species H⁻(t) H2⁺(t) H3⁺(t) CH⁺(t) CH2⁺(t) OH⁺(t) H2O⁺(t) H3O⁺(t) CO⁺(t) HOC⁺(t) O⁻(t) C⁻(t) O2⁺(t) e(t) H⁺(t) H(t) H2(t) He(t) He⁺(t) C(t) C⁺(t) O(t) O⁺(t) OH(t) H2O(t) CO(t) C2(t) O2(t) H2O⁺(t) CH(t) CH2(t) CH3⁺(t) M(t)
@parameters T n_H Av k1 k2 k6 k12 k14 k15 k17 k19 k20 k21 k29 k38 k42 k47 k52 k146 k149 k154 k157 k158 R188 R189 R190 R191 R192 R193 R194 R195


print("Checkpoint 4\n")
# %% Define the network
reaction_equations = [
    (@reaction n_H * k1, H + e --> H⁻),    #k1
    (@reaction n_H * k2, H⁻ + H --> H2 + e), #k2
    (@reaction n_H * 10^(-19.38 - 1.523 * log10(T) + 1.118 * (log10(T))^2 - 0.1269 * (log10(T))^3) , H + H⁺ --> H2⁺ ), 
    (@reaction n_H * 6.4e-10 , H + H2⁺ --> H2 + H⁺), 
    (@reaction n_H * 2.4e-6 * (T)^(-1/2) * (1 + T/20000) , H⁻ + H⁺--> H + H), 
    (@reaction n_H * k6, H2⁺ + e --> H + H), 
    (@reaction n_H * (-3.3232183e-7 + 3.3735382e-7 * log(T) - 1.4491368e-7 * (log(T))^2 + 3.4172805e-8 * (log(T))^3 - 4.781372e-9 * (log(T))^4 + 3.9731542e-10 * (log(T))^5 - 1.8171411e-11 * (log(T))^6 + 3.5311932e-13 * (log(T))^7) * exp(-21237.15/T) , H2 + H⁺ --> H2⁺ + H), 
    (@reaction n_H * 3.73e-9 * (T)^(0.1121) * exp(-99430/T) , H2 + e --> H + H + e), 
    (@reaction n_H * 6.67e-12 * (T)^(1/2) * exp(-(1 + 63590/T)) , H2 + H --> H + H + H), 
    (@reaction n_H * (5.996e-30 * (T)^(4.1881)) / ((1 + 6.761e-6 * T)^(5.6881)) * exp(-54657.4/T) , H2 + H2 --> H2 + H + H), 
    (@reaction n_H * exp(-3.271396786e1 + 1.3536556e1 * log(Te) - 5.73932875 * (log(Te))^2 + 1.56315498 * (log(Te))^3 - 2.877056e-1 * (log(Te))^4 + 3.48255977e-2 * (log(Te))^5 - 2.63197617e-3 * (log(Te))^6 + 1.11954395e-4 * (log(Te))^7 - 2.03914985e-6 * (log(Te))^8) , H + e --> H⁺ + e + e), 
    (@reaction n_H * k12, H⁺ + e --> H ),
    (@reaction n_H * exp(-1.801849334e1 + 2.3608522 * log(Te) - 2.827443e-1 * (log(Te))^2 + 1.62331664e-2 * (log(Te))^3 - 3.36501203e-2 * (log(Te))^4 + 1.17832978e-2 * (log(Te))^5 - 1.6561947e-3 * (log(Te))^6 + 1.0682752e-4 * (log(Te))^7 - 2.63128581e-6 * (log(Te))^8) , H⁻ + e --> H + e + e), 
    (@reaction n_H * k14, H⁻ + H --> H + H + e), 
    (@reaction n_H * k15, H⁻ + H⁺--> H2⁺ + e),
    (@reaction n_H * exp(-4.409864886e1 + 2.391596563e1 * log(Te) - 1.07532302e1 * (log(Te))^2 + 3.05803875 * (log(Te))^3 - 5.6851189e-1 * (log(Te))^4 + 6.79539123e-2 * (log(Te))^5 - 5.0090561e-3 * (log(Te))^6 + 2.06723616e-4 * (log(Te))^7 - 3.64916141e-6 * (log(Te))^8) , He + e --> He⁺ + e + e), 
    (@reaction n_H * k17, He⁺ + e --> He ), 
    (@reaction n_H * 1.25e-15 * (T/300)^(0.25) , He⁺ + H --> He + H⁺), 
    (@reaction n_H * k19, He + H⁺--> He⁺ + H), 

    
    #(@reaction n_H * 4.67e-12 * (T/300)^(-0.6), C⁺ + e --> C ), # Glover rate for C⁺ + e --> C 
    (@reaction n_H * 2.36e-12 * (T/300)^(-0.29) * exp(-(-17.6/T)), C⁺ + e --> C ), # Umist rate for C⁺ + e --> C 
    #(@reaction n_H * 1.4e-10 * (T)^(-0.61), C⁺ + e --> C ), # Nelson rate 1 for C⁺ + e --> C

    (@reaction 2.3e-17, C --> C⁺ + e ), # Umist rate for C --> C⁺ + e
    #(@reaction 3e-10 * 2 * exp(-3 * Av), C --> C⁺ + e ), # Nelson rate 2 for C --> C⁺ + e
    

    
    (@reaction n_H * k21, O⁺ + e --> O ), 
    (@reaction n_H * 6.85e-8 * (0.193 + (11.26/Te))^(-1) * (11.26/Te)^(0.25) * exp(-11.26/Te) , C + e --> C⁺ + e + e), 
    (@reaction n_H * 3.59e-8 * (0.073 + (13.6/Te))^(-1) * (13.6/Te)^(0.34) * exp(-13.6/Te) , O + e --> O⁺ + e + e), 
    (@reaction n_H * (4.99e-11 * T^(0.405) + 7.54e-10 * T^(-0.458)) , O⁺ + H --> O + H⁺), 
    (@reaction n_H * (1.08e-11 * T^(0.517) + 4e-10 * T^(0.00669)) * exp(-227/T) , O + H⁺--> O⁺ + H), 
    (@reaction n_H * (4.991e-15 * (T/10000)^(0.3794) * exp(-T/1121000) + 2.78e-15 * (T/10000)^(-0.2163)) * exp(T/815800) , O + He⁺ --> O⁺ + He), 
    (@reaction n_H * 3.9e-16 * T^(0.213) , C + H⁺ --> C⁺ + H), 
    (@reaction n_H * 6.08e-14 * (T/10000)^(1.96) * exp(-170000/T) , C⁺ + H --> C + H⁺), 
    (@reaction n_H * k29, C + He⁺ --> C⁺ + He), 
    (@reaction n_H * 10^(-27.029 + 3.801 * log10(T) - 29487/T) , H2 + He --> H + H + He), 
    (@reaction n_H * 6e-9 * exp(-50900/T) , OH + H --> O + H + H), 
    (@reaction n_H * 3.8e-10 , HOC⁺ + H2 --> HCO⁺ + H2), 
    (@reaction n_H * 4e-10 , HOC⁺ + CO --> HCO⁺ + CO), 
    (@reaction n_H * 6.64e-10 * exp(-11700/T) , C + H2 --> CH + H), 
    (@reaction n_H * 1.31e-10 * exp(-80/T) , CH + H --> C + H2), 
    (@reaction n_H * 5.46e-10 * exp(-1943/T) , CH + H2 --> CH2 + H), 
    (@reaction n_H * 6.59e-11 , CH + C --> C2 + H), 
    (@reaction n_H * k38, CH+O --> CO+H), # Nelson Match Neutral-Neutral reaction 1, Nelson has 2e-10
    (@reaction n_H * 6.64e-11 , CH2 + H --> CH + H2), 
    (@reaction n_H * 1.33e-10 , CH2 + O --> CO + H + H), 
    (@reaction n_H * 8e-11 , CH2 + O --> CO + H2), 
    (@reaction n_H * k42, C2 + O --> CO + C), 
    (@reaction n_H * 3.14e-13 * (T/300)^(2.7) * exp(-3150/T) , O + H2 --> OH + H), 
    (@reaction n_H * 6.99e-14 * (T/300)^(2.8) * exp(-1950/T) , OH + H --> O + H2), 
    (@reaction n_H * 2.05e-12 * (T/300)^(1.52) * exp(-1736/T) , OH + H2 --> H2O + H), 
    (@reaction n_H * 1e-10 , OH + C --> CO + H), # Nelson Match Neutral-Neutral reaction 2, Nelson has 5.8e-12 * T^(0.5)
    (@reaction n_H * k47, OH + O --> O2 + H), 
    (@reaction n_H * 1.65e-12 * (T/300)^(1.14) * exp(-50/T) , OH + OH --> H2O + H), 
    (@reaction n_H * 1.59e-11 * (T/300)^(1.2) * exp(-9610/T) , H2O + H --> H2 + OH), 
    (@reaction n_H * 2.61e-10 * exp(-8156/T) , O2 + H --> OH + O), 
    (@reaction n_H * 3.16e-10 * exp(-21890/T) , O2 + H2 --> OH + OH), 
    (@reaction n_H * k52, O2 + C --> CO + O), 
    (@reaction n_H * 1.1e-10 * (T/300)^(0.5) * exp(-77700/T) , CO + H --> C + OH), 
    (@reaction n_H * 2.24e-9 * (T/300)^(0.042) * exp(-T/46600) , H2⁺ + H2 --> H3⁺ + H), 
    (@reaction n_H * 7.7e-9 * exp(-17560/T) , H3⁺ + H --> H2⁺ + H2), 
    (@reaction n_H * 2.4e-9 , C + H2⁺ --> CH⁺ + H), 
    (@reaction n_H * 2e-9 , C + H3⁺ --> CH⁺ + H2), # Nelson Match Ion-Molecule reaction 1
    (@reaction n_H * 1e-10 * exp(-4640/T) , C⁺ + H2 --> CH⁺ + H), # Nelson Match Ion-Molecule reaction 6, Nelson has 4e-16
    (@reaction n_H * 7.5e-10 , CH⁺ + H --> C⁺ + H2), 
    (@reaction n_H * 1.2e-9 , CH⁺ + H2 --> CH2⁺ + H), 
    (@reaction n_H * 3.5e-10 , CH⁺ + O --> CO⁺ + H), 
    (@reaction n_H * 1.4e-9 , CH2 + H⁺--> CH⁺ + H2), 
    (@reaction n_H * 1e-9 * exp(-7080/T) , CH2⁺ + H --> CH⁺ + H2), 
    (@reaction n_H * 1.6e-9 , CH2⁺ + H2 --> CH3⁺ + H), 
    (@reaction n_H * 7.5e-9 , CH2⁺ + O --> HCO⁺ + H), 
    (@reaction n_H * 7e-10 * exp(-10560/T) , CH3⁺ + H --> CH2⁺ + H2), 
    (@reaction n_H * 4e-10 , CH3⁺ + O --> HCO⁺ + H2), 
    (@reaction n_H * 4.8e-10 , C2 + O⁺ --> CO⁺ + C), 
    (@reaction n_H * 1.7e-9 , O⁺ + H2 --> OH⁺ + H), 
    (@reaction n_H * 1.5e-9 , O + H2⁺ --> OH⁺ + H), 
    (@reaction n_H * 8.4e-10 , O + H3⁺ --> OH⁺ + H2), # Nelson Match Ion-Molecule reaction 2 Nelson has 8e-10
    (@reaction n_H * 1.3e-9 , OH + H3⁺ --> H2O⁺ + H2), 
    (@reaction n_H * 7.7e-10 , OH + C⁺ --> CO⁺ + H), 
    (@reaction n_H * 1.01e-9 , OH⁺ + H2 --> H2O⁺ + H), 
    (@reaction n_H * 6.4e-10 , H2O⁺ + H2 --> H3O⁺ + H), 
    (@reaction n_H * 5.9e-9 , H2O + H3⁺ --> H3O⁺ + H2), 
    (@reaction n_H * 9e-10 , H2O + C⁺ --> HCO⁺ + H), 
    (@reaction n_H * 1.8e-9 , H2O + C⁺ --> HOC⁺ + H), 
    (@reaction n_H * 1e-11 , H3O⁺ + C --> HCO⁺ + H2), 
    (@reaction n_H * 3.8e-10 , O2 + C⁺ --> CO⁺ + O), 
    (@reaction n_H * 6.2e-10 , O2 + C⁺ --> CO + O⁺), 
    (@reaction n_H * 9.1e-10 , O2 + CH2⁺ --> HCO⁺ + OH), 
    (@reaction n_H * 5.2e-11 , O2⁺ + C --> CO⁺ + O), 
    (@reaction n_H * 2.7e-11 , CO + H3⁺ --> HOC⁺ + H2), 
    (@reaction n_H * 1.7e-9 , CO + H3⁺ --> HCO⁺ + H2), # Nelson Match Ion-Molecule reaction 3
    (@reaction n_H * 1.1e-9 , HCO⁺ + C --> CO + CH⁺), 
    (@reaction n_H * 2.5e-9 , HCO⁺ + H2O --> CO + H3O⁺), 
    (@reaction n_H * 7.2e-15 , H2 + He⁺  -->  He + H2⁺), # Nelson Match Ion-Molecule reaction 4, Nelson has 7e-15
    (@reaction n_H * 3.7e-14 * exp(-35/T) , H2 + He⁺ --> He + H + H⁺), # Nelson Match Ion-Molecule reaction 4, Nelson has 7e-15
    (@reaction n_H * 1.9e-9 , CH + H⁺ -->  CH⁺ + H), 
    (@reaction n_H * 1.4e-9 , CH2 + H⁺ -->  CH2⁺ + H), 
    (@reaction n_H * 7.5e-10 , CH2 + He⁺ -->  C⁺ + He + H2), 
    (@reaction n_H * 1.6e-9 , C2 + He⁺ -->  C⁺ + C + He), 
    (@reaction n_H * 2.1e-9 , OH + H⁺ -->  OH⁺ + H), 
    (@reaction n_H * 1.1e-9 , OH + He⁺ -->  O⁺ + He + H), 
    (@reaction n_H * 6.9e-9 , H2O + H⁺ -->  H2O⁺ + H), 
    (@reaction n_H * 2.04e-10 , H2O + He⁺ -->  OH + He + H⁺), 
    (@reaction n_H * 2.86e-10 , H2O + He⁺ -->  OH⁺ + He + H), 
    (@reaction n_H * 6.05e-11 , H2O + He⁺ -->  H2O⁺ + He), 
    (@reaction n_H * 2e-9 , O2 + H⁺ -->  O2⁺ + H), 
    (@reaction n_H * 3.3e-11 , O2 + He⁺ -->  O2⁺ + He), 
    (@reaction n_H * 1.1e-9 , O2 + He⁺ -->  O⁺ + O + He), 
    (@reaction n_H * 5.2e-11 , O2⁺ + C  -->  O2 + C⁺), 
    (@reaction n_H * 1.4e-9*(T/300)^(-0.5) , CO + He⁺ -->  C⁺ + O + He), # Nelson Match Ion-Molecule reaction 5, Nelson has 1.6e-9
    (@reaction n_H * 1.4e-16 * (T/300)^(-0.5) , CO + He⁺ -->  C + O⁺ + He), 
    (@reaction n_H * 7.5e-10 , CO⁺ + H  -->  CO + H⁺), 
    (@reaction n_H * 2.3e-7 * (T/300)^(-0.5) , C⁻ + H⁺ -->  C + H), 
    (@reaction n_H * 2.3e-7 * (T/300)^(-0.5) , O⁻ + H⁺ -->  O + H), 
    (@reaction n_H * 2.32e-7 * (T/300)^(-0.52) * exp(T/22400) , He⁺ + H⁻  -->  He + H), 
    (@reaction n_H * 2.34e-8 * (T/300)^(-0.52) , H3⁺ + e --> H2 + H), 
    (@reaction n_H * 4.36e-8 * (T/300)^(-0.52) , H3⁺ + e --> H + H + H), 
    (@reaction n_H * 7e-8 * (T/300)^(-0.5) , CH⁺ + e --> C + H), 
    (@reaction n_H * 1.6e-7 * (T/300)^(-0.6) , CH2⁺ + e --> CH + H), 
    (@reaction n_H * 4.03e-7 * (T/300)^(-0.6) , CH2⁺ + e --> C + H + H), 
    (@reaction n_H * 7.68e-8 * (T/300)^(-0.6) , CH2⁺ + e --> C + H2), 
    (@reaction n_H * 7.75e-8 * (T/300)^(-0.5) , CH3⁺ + e --> CH2 + H), 
    (@reaction n_H * 1.95e-7 * (T/300)^(-0.5) , CH3⁺ + e --> CH + H2), 
    (@reaction n_H * 2e-7 * (T/300)^(-0.4) , CH3⁺ + e --> CH + H + H), 
    (@reaction n_H * 6.3e-9 * (T/300)^(-0.48) , OH⁺ + e --> O + H), 
    (@reaction n_H * 3.05e-7 * (T/300)^(-0.5) , H2O⁺ + e --> O + H + H), 
    (@reaction n_H * 3.9e-8 * (T/300)^(-0.5) , H2O⁺ + e --> O+H2), 
    (@reaction n_H * 8.6e-8 * (T/300)^(-0.5) , H2O⁺ + e --> OH + H), 
    (@reaction n_H * 1.08e-7 * (T/300)^(-0.5) , H3O⁺ + e --> H + H2O), 
    (@reaction n_H * 6.02e-8 * (T/300)^(-0.5) , H3O⁺ + e --> OH + H2), 
    (@reaction n_H * 2.58e-7 * (T/300)^(-0.5) , H3O⁺ + e --> OH + H + H), 
    (@reaction n_H * 5.6e-9 * (T/300)^(-0.5) , H3O⁺ + e --> O + H + H2), 
    (@reaction n_H * 1.95e-7 * (T/300)^(-0.7) , O2⁺ + e --> O + O), 
    (@reaction n_H * 2.75e-7 * (T/300)^(-0.55) , CO⁺ + e --> C + O), 
    (@reaction n_H * 2.76e-7 * (T/300)^(-0.64) , HCO⁺ + e --> CO + H), 
    (@reaction n_H * 2.4e-8 * (T/300)^(-0.64) , HCO⁺ + e --> OH + C), 
    (@reaction n_H * 1.1e-7 * (T/300)^(-1) , HOC⁺ + e --> CO + H), 
    (@reaction n_H * 1e-9 , H⁻ + C --> CH + e), 
    (@reaction n_H * 1e-9 , H⁻ + O --> OH + e), 
    (@reaction n_H * 1e-10 , H⁻ + OH --> H2O + e), 
    (@reaction n_H * 5e-10 , C⁻ + H --> CH + e), 
    (@reaction n_H * 1e-13 , C⁻ + H2 --> CH2 + e), 
    (@reaction n_H * 5e-10 , C⁻ + O -->  CO + e), 
    (@reaction n_H * 5e-10 , O⁻ + H  -->  OH + e), 
    (@reaction n_H * 7e-10 , O⁻ + H2 -->  H2O + e), 
    (@reaction n_H * 5e-10 , O⁻ + C --> CO + e), 
    (@reaction n_H * 1e-16 , H2 + H⁺ --> H3⁺), 
    (@reaction n_H * 2.25e-15 , C + e --> C⁻ ), 
    (@reaction n_H * 1e-17 , C + H --> CH ), 
    (@reaction n_H * 1e-17 , C + H2 --> CH2 ), 
    (@reaction n_H * 4.36e-18 * (T/300)^(0.35) * exp(-161.3/T) , C + C --> C2 ), 
    (@reaction n_H * k146, C + O  --> CO ), 
    (@reaction n_H * 4.46e-16 * T^(-0.5) * exp(-4.93/(T^(2/3))) , C⁺ + H --> CH⁺ ), 
    (@reaction n_H * 4e-16 *(T/300)^(-0.2) , C⁺ + H2 --> CH2⁺ ), 
    (@reaction n_H * k149, C⁺ + O --> CO⁺ ), 
    (@reaction n_H * 1.5e-15 , O + e --> O⁻ ), 
    (@reaction n_H * 9.9e-19 * (T/300)^(-0.38) , O+H --> OH ), 
    (@reaction n_H * 4.9e-20 * (T/300)^1.58 , O + O --> O2 ), 
    (@reaction n_H * 5.26e-18 *(T/300)^(-5.22) * exp(-90/T) , OH + H --> H2O ), 
    (@reaction n_H^2 * k154, H + H + H --> H2 + H), 
    (@reaction n_H^2 * 2.8e-31 * (T)^-0.6 , H + H + H2 --> H2 + H2), 
    (@reaction n_H^2 * 6.9e-32 * (T)^(-0.4) , H + H + He --> H2 + He), 
    (@reaction n_H^2 * k157, C + C + M --> C2 + M), 
    (@reaction n_H^2 * k158, C + O + M --> CO + M), 
    (@reaction n_H^2 * 100 * cr_ion_rate * 960, C⁺ + O + M --> CO⁺ + M), # FIX: these rely on k210, which doesn't exist?
    (@reaction n_H^2 * 100 * cr_ion_rate * 960, C + O⁺ + M --> CO⁺ + M), # FIX: these rely on k210, which doesn't exist?
    (@reaction n_H^2 * 4.33e-32 * (T/300)^(-1) , O + H + M --> OH + M), 
    (@reaction n_H^2 * 2.56e-31 * (T/300)^(-2) , OH + H + M --> H2O + M), 
    (@reaction n_H^2 * 9.2e-34 * (T/300)^(-1) , O + O + M --> O2 + M), 
    (@reaction n_H * 2e-11 * (T/300)^(0.44) , O + CH -->  HCO⁺ + e), 
    (@reaction n_H * 3e-18 * (T)^(0.5) * (1 + 1e4 * exp(-600/Td))^(-1) * (1 + 0.04*(T + Td)^(0.5) + 0.002*T + 8e-6 * T^2)^(-1) , H+H --> H2), # FIX H is H(s) here?
    


    # Photochemical Reactions (R166) start below
    (@reaction 7.1e-7 * exp(-0.5 * Av), H⁻ --> H + e), 
    (@reaction 1.1e-9 * exp(-1.9 * Av), H2⁺ --> H + H⁺), 
    (@reaction 5.6e-11 * exp(-2.55*Av + 0.0165*(Av)^2), H2   -->  H+H), 
    (@reaction 4.9e-13 * exp(-1.8 * Av), H3⁺ --> H2 + H⁺),
    (@reaction 4.9e-13 * exp(-2.3 * Av), H3⁺ --> H2⁺ + H), 
    (@reaction 3.1e-10 * exp(-3.0 * Av), C -->  C⁺ + e), #171
    (@reaction 2.4e-7 * exp(-0.9 * Av), C⁻ --> C + e), 
    (@reaction 8.7e-10 * exp(-1.2 * Av), CH --> C + H), 
    (@reaction 7.7e-10 * exp(-2.8 * Av), CH --> CH⁺ + e), 
    (@reaction 2.6e-10 * exp(-2.5 * Av), CH⁺ --> C + H⁺), 
    (@reaction 7.1e-10 * exp(-1.7 * Av), CH2 --> CH + H), 
    (@reaction 5.9e-10 * exp(-2.3 * Av), CH2 --> CH2⁺ + e), 
    (@reaction 4.6e-10 * exp(-1.7 * Av), CH2⁺ --> CH⁺ + H), 
    (@reaction 1.0e-9 * exp(-1.7 * Av), CH3⁺ --> CH2⁺ + H), 
    (@reaction 1.0e-9 * exp(-1.7 * Av), CH3⁺ --> CH⁺ + H2), 
    (@reaction 1.5e-10 * exp(-2.1 * Av), C2 --> C + C), 
    (@reaction 2.4e-7 * exp(-0.5 * Av), O⁻ --> O + e),
    (@reaction 3.7e-10 * exp(-1.7 * Av), OH --> O + H), 
    (@reaction 1.6e-12 * exp(-3.1 * Av), OH --> OH⁺ + e), 
    (@reaction 1.0e-12 * exp(-1.8 * Av), OH⁺  --> O + H⁺), 
    (@reaction 6.0e-10 * exp(-1.7 * Av), H2O --> OH + H), 
    (@reaction 3.2e-11 * exp(-3.9 * Av), H2O --> H2O⁺ + e), 
    (@reaction R188, H2O⁺ --> H2⁺ + O), 
    (@reaction R189, H2O⁺ --> H⁺ + OH), 
    (@reaction R190, H2O⁺ --> O⁺ + H2), 
    (@reaction R191, H2O⁺ --> OH⁺ + H), 
    (@reaction R192, H3O⁺ --> H⁺ + H2O), 
    (@reaction R193, H3O⁺ --> H2⁺ + OH), 
    (@reaction R194, H3O⁺ --> H2O⁺ + H), 
    (@reaction R195, H3O⁺ --> OH⁺ + H2), 
    (@reaction 5.6e-11 * exp(-3.7 * Av), O2 --> O2⁺ + e), 
    (@reaction 7.0e-10 * exp(-1.8 * Av), O2 --> O + O), 
    (@reaction 2e-10 * exp(-3.5 * Av), CO   --> C + O), #from UMIST


    # Cosmic Ray Reactions (R199) start below 
    (@reaction cr_ion_rate * 1 , H --> H⁺ + e), 
    (@reaction cr_ion_rate * 1.1333333 , He --> He⁺ + e), # cr_rate = 6e-18, Nelson match CR Ionization reaction 2
    (@reaction cr_ion_rate * 0.037 , H2 --> H⁺ + H + e), 
    (@reaction cr_ion_rate * 0.22 , H2 --> H + H), 
    (@reaction cr_ion_rate * 6.5e-4 , H2 --> H⁺ + H⁻), 
    (@reaction cr_ion_rate * 2 , H2 --> H2⁺ + e),  # cr_rate = 6e-18, Nelson match CR Ionization reaction 1
    (@reaction cr_ion_rate * 3.8 , C --> C⁺ + e), 
    (@reaction cr_ion_rate * 5.7 , O --> O⁺ + e), 
    (@reaction cr_ion_rate * 6.5 , CO --> CO⁺ + e), 
    (@reaction cr_ion_rate * 2800 , C --> C⁺ + e), 
    (@reaction cr_ion_rate * 4000 , CH --> C + H), 
    (@reaction cr_ion_rate * 960 , CH⁺ --> C⁺ + H), 
    (@reaction cr_ion_rate * 2700 , CH2 --> CH2⁺ + e), 
    (@reaction cr_ion_rate * 2700 , CH2 --> CH + H), 
    (@reaction cr_ion_rate * 1300 , C2 --> C + C), 
    (@reaction cr_ion_rate * 2800 , OH --> O + H), 
    (@reaction cr_ion_rate * 5300 , H2O --> OH + H), 
    (@reaction cr_ion_rate * 4100 , O2 --> O + O), 
    (@reaction cr_ion_rate * 640  , O2 --> O2⁺ + e), 
    #(@reaction cr_ion_rate * 2e-10 * exp(-3.5 * Av), CO   --> C + O) # this is a cosmic ray reaction, NOT photoreaction
    #(@reaction 10e-10 * (1) * (1.7) * exp(-3 *  Av) , CO --> C + O) # I copied the rate for this one from nelson
]

print("Checkpoint 5 \n")

    # %% Turn the network into an ODE system
@named system = ReactionSystem(reaction_equations, t)
odesys = convert(ODESystem, complete(system))
sys = convert(ODESystem, complete(system))
ssys = structural_simplify(sys)
prob = ODEProblem(ssys, u0, tspan, params)
#sol = solve(prob, lsoda(), reltol=1.49012e-8, abstol=1.49012e-8, saveat=1e10)
#sol = solve(prob, lsoda(), saveat=1e3)
sol = solve(prob, reltol=1.49012e-6, abstol=1.49012e-6, Rodas4())


# C and C+
plot(sol, idxs = (0,10), lw = 3, lc = "blue", title = "Glover with Glover rates")
plot!(sol, idxs = (0,9), lw = 3, lc = "orange", title = "Glover with Glover rates")

#=

# CO
plot(sol, idxs = (0,16), lw = 3, lc = "green", title = "Glover")

# O 
plot(sol, idxs = (0,12), lw = 3, lc = "blue", title = "Glover")

# He+
plot(sol, idxs = (0,8), lw = 3, lc = "light pink", title = "Glover")

# CH and CH2 = CHx
plot(sol, idxs = (0,17), lw = 3, lc = "blue", title = "Glover")
plot!(sol, idxs = (0,18), lw = 3, lc = "light blue", title = "Glover")

# OH, OH+, H2O, H2O+, and O2 = OHx
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





