Updated the nelson ODE function, (made it nicer to look at) and added parameters

    Hello! My name is Nina De La Torre, and I wrote this ODE function! It was such a pleasant surprise to see it on the Julia documentation website, but the version that was used was unfinished, and I wanted to fix it up a bit. Anyways, here is the same function but a bit nicer to look at, with helpful comments, parameters, and a cute name. I am planning to request updates the other nelson files as well. I hope this is helpful!

    Kindly, 
    Nina









    … parameters

## Checklist

- [x] Appropriate tests were added
- [x] Any code changes were done in a way that does not break public API
- [x] All documentation related to code changes were updated
- [x] The new code follows the
  [contributor guidelines](https://github.com/SciML/.github/blob/master/CONTRIBUTING.md), in particular the [SciML Style Guide](https://github.com/SciML/SciMLStyle) and
  [COLPRAC](https://github.com/SciML/COLPRAC).
- [x] Any new documentation only uses public API
  
## Additional context

Here are the Nelson ODE updates! Feel free to take or leave any portion of these changes. I didn't find anywhere to add figures, but the plot that this code produces just shows a funky comparison of different algorithms. Back when I was searching for accurate, quick algorithms to use
Thanks!






#using CUDA, LinearAlgebra

 @time begin
    # %% Set The timespan, paramters, and initial conditions
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs

params = Dict(
    :T => 10,  
    :n_H => 611 )

u0 = [  
    1, # 1: H_m, H-
    1, # 2: H2_p
    1, # 3: H3_p
    1, # 4: CH_p
    1, # 5: CH2_p
    1, # 6: OH_p
    1, # 7: H2O_p
    1, # 8: H3O_p
    1, # 9: CO_p
    1, # 10: HOC_p
    1, # 11: O_m
    1, # 12: C_m
    1, # 13: O2_p
    1, # 14: e
    1, # 15: H_p
    1, # 16: H
    1, # 17: H2
    1, # 18: He
    1, # 19: He_p
    1, # 20: C
    1, # 21: C_p
    1, # 22: O
    1, # 23: O_p
    1, # 24: OH
    1, # 25: H2O
    1, # 26: CO
    1, # 27: C2
    1, # 28: O2
    1, # 29: H2O_p
    1, # 30: CH
    1, # 31: CH2
    1, # 32: CH3_p
    1, # 33: M
    ]


    # %% Network admin things
allvars = @strdict params tspan u0
@variables t 
@species H_m(t) H2_p(t) H3_p(t) CH_p(t) CH2_p(t) OH_p(t) H2O_p(t) H3O_p(t) CO_p(t) HOC_p(t) O_m(t) C_m(t) O2_p(t) e(t) H_p(t) H(t) H2(t) He(t) He_p(t) C(t) C_p(t) O(t) O_p(t) OH(t) H2O(t) CO(t) C2(t) O2(t) HCO_p(t) CH(t) CH2(t) CH3_p(t) M(t)
@parameters T n_H
#CHECK I am replacing ln with log

    # %% Define the network
reaction_equations = [
    (@reaction 10^(-17.845 + 0.762 * log(T) + 0.1523 * (log(T))^2 - 0.03274 * (log(T))^3) , H+e --> H_m), 
    (@reaction 0.0000000015 , H_m+H --> H2+e), 
    (@reaction 10^(-19.38 - 1.523 * log(T) + 1.118 * (log(T))^2 - 0.1269 * (log(T))^3) , H+ H_p--> H2_p ), 
    (@reaction 0.00000000064 , H+H2_p --> H2 + H_p), 
    (@reaction 0.0000024 * (T)^(-1/2) * (1 + T/20000) , H_m + H_p--> H+H), 
    (@reaction 0.00000001 , H2+ + e --> H+H), 
    (@reaction (-0.00000033232183 + 0.00000033735382 * log(T) - 0.00000014491368 * (log(T))^2 + 0.000000034172805 * (log(T))^3 - 0.000000004781372 * (log(T))^4 + 0.00000000039731542 * (log(T))^5 - 0.000000000018171411 * (log(T))^6 + 35311932000000 * (log(T))^7) * exp(-21237.15/T) , H2+H_p--> H2_p+H), 
    (@reaction 0.00000000373 * (T)^(0.1121) * exp(-99430/T) , H2+e --> H+H+e), 
    (@reaction 0.00000000000667 * (T)^(1/2) * exp(-(1 + 63590/T)) , H2+H --> H+H+H), 
    (@reaction (5.996E+30 * (T)^(4.1881)) / ((1 + 0.000006761 * T)^(5.6881)) * exp(-54657.4/T) , H2+H2 --> H2+H+H), 
    (@reaction exp(-32.71396786 + 13.536556 * log(T) - 5.73932875 * (log(T))^2 + 1.56315498 * (log(T))^3 - 0.2877056 * (log(T))^4 + 0.0348255977 * (log(T))^5 - 0.00263197617 * (log(T))^6 + 0.000111954395 * (log(T))^7 - 0.00000203914985 * (log(T))^8) , H+e --> H_p + e+e), 
    (@reaction 0.0000000000001269 * ((315614/T)^(1.503)) * (1 + (604625/T)^(0.47))^(-1.923) , H_p + e --> H ),
    (@reaction exp(-18.01849334 + 2.3608522 * log(T) - 0.2827443 * (log(T))^2 + 0.0162331664 * (log(T))^3 - 0.0336501203 * (log(T))^4 + 0.0117832978 * (log(T))^5 - 0.0016561947 * (log(T))^6 + 0.00010682752 * (log(T))^7 - 0.00000263128581 * (log(T))^8) , H_m + e --> H + e + e), 
    (@reaction 0.0000000025634 * (T)^(1.78186) , H_m + H --> H + H + e), 
    (@reaction 0.0000000069 * T^(-0.35) , H_m + H_p--> H2_p+e),
    (@reaction exp(-44.09864886 + 23.91596563 * log(T) - 10.7532302 * (log(T))^2 + 3.05803875 * (log(T))^3 - 0.56851189 * (log(T))^4 + 0.0679539123 * (log(T))^5 - 0.0050090561 * (log(T))^6 + 0.000206723616 * (log(T))^7 - 0.00000364916141 * (log(T))^8) , He+e --> He_p + e+e), 
    (@reaction 0.0000000001 * T^(-0.5) * (12.72 - 1.615 * log(T) - 0.3162 * (log(T))^2 + 0.0493 * (log(T))^3) , He_p + e --> He ), 
    (@reaction 0.00000000000000125 * (T/300)^(0.25) , He_p + H --> He+H_p), 
    (@reaction 0.00000000126 * T^(-0.75) * exp(-127500/T) , He+H_p--> He_p + H), 
    (@reaction 0.00000000000467 * (T/300)^(-0.6) , C_p + e --> C ), 
    (@reaction 0.00000000013 * T^(-0.64) , O_p + e --> O ), 
    (@reaction 0.0000000685 * (0.193 + (11.26/T))^(-1) * (11.26/T)^(0.25) * exp(-(11.26/T)) , C+e --> C_p + e+e), 
    (@reaction 0.0000000359 * (0.073 + (13.6/T))^(-1) * (13.6/T)^(0.34) * exp(-13.6/T) , O+e --> O_p + e+e), 
    (@reaction 0.0000000000499 * T^(0.405) + 0.000000000754 * T^(-0.458) , O_p + H --> O+H_p), 
    (@reaction (0.0000000000108 * T^(0.517) + 0.0000000004 * T^(0.00669)) * exp(-227/T) , O+H_p--> O_p + H), 
    (@reaction 0.000000000000004991 * (T/10000)^(0.3794) * exp(-T/1121000) + 0.00000000000000278 * (T/10000)^(-0.2163) * exp(T/815800) , O+He_p --> O_p + He), 
    (@reaction 0.00000000000000039 * T^(0.213) , C+H_p--> C_p + H), 
    (@reaction 0.0000000000000608 * (T/10000)^(1.96) * exp(-170000/T) , C_p + H --> C+H_p), 
    (@reaction 0.0000000000000000858 * (T)^(0.757) , C+He_p --> C_p + He), 
    (@reaction 10^(-27.029 + 3.801 * log(T) - 29487/T) , H2+He --> H+H+He), 
    (@reaction 0.000000006 * exp(-50900/T) , OH+H --> O+H+H), 
    (@reaction 0.00000000038 , HOC_p + H2 --> HCO_p + H2), 
    (@reaction 0.0000000004 , HOC_p + CO --> HCO_p + CO), 
    (@reaction 0.000000000664 * exp(-11700/T) , C+H2 --> CH+H), 
    (@reaction 0.000000000131 * exp(-80/T) , CH+H --> C+H2), 
    (@reaction 0.000000000546 * exp(-1943/T) , CH+H2 --> CH2+H), 
    (@reaction 0.0000000000659 , CH+C --> C2+H), 
    (@reaction 0.000000000066 , CH+O --> CO+H), 
    (@reaction 0.0000000000664 , CH2+H --> CH+H2), 
    (@reaction 0.000000000133 , CH2+O --> CO+H+H), 
    (@reaction 0.00000000008 , CH2+O --> CO+H2), 
    (@reaction 0.00000000005 * (T/300)^(0.5) , C2+O --> CO+C), 
    (@reaction 0.000000000000314 * (T/300)^(2.7) * exp(-3150/T) , O+H2 --> OH+H), 
    (@reaction 0.0000000000000699 * (T/300)^(2.8) * exp(-1950/T) , OH+H --> O+H2), 
    (@reaction 0.00000000000205 * (T/300)^(1.52) * exp(-1736/T) , OH+H2 --> H2O+H), 
    (@reaction 0.0000000001 , OH+C --> CO+H), 
    (@reaction 0.000000000035 , OH+O --> O2+H), 
    (@reaction 0.00000000000165 * (T/300)^(1.14) * exp(-50/T) , OH+OH --> H2O+H), 
    (@reaction 0.0000000000159 * (T/300)^(1.2) * exp(-9610/T) , H2O+H --> H2+OH), 
    (@reaction 0.000000000261 * exp(-8156/T) , O2+H --> OH+O), 
    (@reaction 0.000000000316 * exp(-21890/T) , O2+H2 --> OH+OH), 
    (@reaction 0.000000000047 * (T/300)^(-0.34) , O2+C --> CO+O), 
    (@reaction 0.00000000011 * (T/300) * exp(-77700/T) , CO+H --> C+OH), 
    (@reaction 0.00000000224 * (T/300)^(0.042) * exp(-T/46600) , H2_p+H2 --> H3_p+H), 
    (@reaction 0.0000000077 * exp(-17560/T) , H3_p+H --> H2_p+H2), 
    (@reaction 0.0000000024 , C+H2_p --> CH_p + H), 
    (@reaction 0.000000002 , C+H3_p --> CH_p + H2), 
    (@reaction 0.0000000001 * exp(-4640/T) , C_p + H2 --> CH_p + H), 
    (@reaction 0.00000000075 , CH_p + H --> C_p + H2), 
    (@reaction 0.0000000012 , CH_p + H2 --> CH2_p+H), 
    (@reaction 0.00000000035 , CH_p + O --> CO_p + H), 
    (@reaction 0.0000000014 , CH2+H_p--> CH_p + H2), 
    (@reaction 0.000000001 * exp(-7080/T) , CH2_p+H --> CH_p + H2), 
    (@reaction 0.0000000016 , CH2_p+H2 --> CH3_p+H), 
    (@reaction 0.00000000075 , CH2_p+O --> HCO_p + H), 
    (@reaction 0.0000000007 * exp(-10560/T) , CH3_p+H --> CH2_p+H2), 
    (@reaction 0.0000000004 , CH3_p+O --> HCO_p + H2), 
    (@reaction 0.00000000048 , C2+O_p --> CO_p + C), 
    (@reaction 0.0000000017 , O_p + H2 --> OH_p + H), 
    (@reaction 0.0000000015 , O+H2_p --> OH_p + H), 
    (@reaction 0.00000000084 , O+H3_p --> OH_p + H2), 
    (@reaction 0.0000000013 , OH+H3_p --> H2O_p + H2), 
    (@reaction 0.00000000077 , OH+C_p --> CO_p + H), 
    (@reaction 0.00000000101 , OH_p + H2 --> H2O_p + H), 
    (@reaction 0.00000000064 , H2O_p + H2 --> H3O_p + H), 
    (@reaction 0.0000000059 , H2O+H3_p --> H3O_p + H2), 
    (@reaction 0.0000000009 , H2O+C_p --> HCO_p + H), 
    (@reaction 0.0000000018 , H2O+C_p --> HOC_p + H), 
    (@reaction 0.00000000001 , H3O_p + C --> HCO_p + H2), 
    (@reaction 0.00000000038 , O2+C_p --> CO_p + O), 
    (@reaction 0.00000000062 , O2+C_p --> CO+O_p), 
    (@reaction 0.00000000091 , O2+CH2_p --> HCO_p + OH), 
    (@reaction 0.000000000052 , O2_p+C --> CO_p + O), 
    (@reaction 0.000000000027 , CO+H3_p --> HOC_p + H2), 
    (@reaction 0.0000000017 , CO+H3_p --> HCO_p + H2), 
    (@reaction 0.0000000011 , HCO_p + C --> CO+CH_p), 
    (@reaction 0.0000000025 , HCO_p + H2O --> CO+H3O_p), 
    (@reaction 0.0000000000000072 , H2+He_p  -->  He+H2_p), 
    (@reaction 0.000000000000037 * exp(-35/T) , H2+He_p  -->  He+H+H_p), 
    (@reaction 0.0000000019 , CH+H_p -->  CH_p + H), 
    (@reaction 0.0000000014 , CH2+H_p -->  CH2_p+H), 
    (@reaction 0.00000000075 , CH2+He_p -->  C_p + He+H2), 
    (@reaction 0.0000000016 , C2+He_p -->  C_p + C+He), 
    (@reaction 0.0000000021 , OH+H_p -->  OH_p + H), 
    (@reaction 0.0000000011 , OH+He_p -->  O_p + He+H), 
    (@reaction 0.0000000069 , H2O+H_p -->  H2O_p + H), 
    (@reaction 0.000000000204 , H2O+He_p -->  OH+He+H_p), 
    (@reaction 0.000000000286 , H2O+He_p -->  OH_p + He+H), 
    (@reaction 0.0000000000605 , H2O+He_p -->  H2O_p + He), 
    (@reaction 0.000000002 , O2+H_p -->  O2_p+H), 
    (@reaction 0.000000000033 , O2+He_p -->  O2_p+He), 
    (@reaction 0.0000000011 , O2+He_p -->  O_p + O+He), 
    (@reaction 0.000000000052 , O2_p+C  -->  O2+C_p), 
    (@reaction 10000-9*(T/300)^(-0.5) , CO+He_p -->  C_p + O+He), 
    (@reaction 0.00000000000000014 * (T/300)^(-0.5) , CO+He_p -->  C+O_p + He), 
    (@reaction 0.00000000075 , CO_p + H  -->  CO+H_p), 
    (@reaction 0.00000023 * (T/300)^(-0.5) , C_m + H_p -->  C + H), 
    (@reaction 0.00000023 * (T/300)^(-0.5) , O_m + H_p -->  O + H), 
    (@reaction 0.000000232 * (T/300)^(-0.52) * exp(T/22400) , He_p + H_m  -->  He+H), 
    (@reaction 0.0000000234 * (T/300)^(-0.52) , H3_p+e  -->  H2+H), 
    (@reaction 0.0000000436 * (T/300)^(-0.52) , H3_p+e  -->  H+H+H), 
    (@reaction 0.00000007 * (T/300)^(-0.5) , CH_p + e  -->  C+H), 
    (@reaction 0.00000016 * (T/300)^(-0.6) , CH2_p+e  -->  CH+H), 
    (@reaction 0.000000403 * (T/300)^(-0.6) , CH2_p+e  -->  C+H+H), 
    (@reaction 0.0000000768 * (T/300)^(-0.6) , CH2_p+e  -->  C+H2), 
    (@reaction 0.0000000775 * (T/300)^(-0.5) , CH3_p+e  -->  CH2+H), 
    (@reaction 0.000000195 * (T/300)^(-0.5) , CH3_p+e  -->  CH+H2), 
    (@reaction 0.0000002 * (T/300)^(-0.4) , CH3_p+e  -->  CH+H+H), 
    (@reaction 0.0000000063 * (T/300)^(-0.48) , OH_p + e  -->  O+H), 
    (@reaction 0.000000305 * (T/300)^(-0.5) , H2O_p + e  -->  O+H+H), 
    (@reaction 0.000000039 * (T/300)^(-0.5) , H2O_p + e  -->  O+H2), 
    (@reaction 0.000000086 * (T/300)^(-0.5) , H2O_p + e  -->  OH+H), 
    (@reaction 0.000000108 * (T/300)^(-0.5) , H3O_p + e  -->  H+H2O), 
    (@reaction 0.0000000602 * (T/300)^(-0.5) , H3O_p + e  -->  OH+H2), 
    (@reaction 0.000000258 * (T/300)^(-0.5) , H3O_p + e  -->  OH+H+H), 
    (@reaction 0.0000000056 * (T/300)^(-0.5) , H3O_p + e  -->  O+H+H2), 
    (@reaction 0.000000195 * (T/300)^(-0.7) , O2_p+e  -->  O+O), 
    (@reaction 0.000000275 * (T/300)^(-0.55) , CO_p + e  -->  C+O), 
    (@reaction 0.000000276 * (T/300)^(-0.64) , HCO_p + e  -->  CO+H), 
    (@reaction 0.000000024 * (T/300)^(-0.64) , HCO_p + e  -->  OH+C), 
    (@reaction 0.00000011 * (T/300)^(-1) , HOC_p + e  -->  CO+H), 
    (@reaction 0.000000001 , H_m+C  -->  CH+e), 
    (@reaction 0.000000001 , H_m+O  -->  OH+e), 
    (@reaction 0.0000000001 , H_m+OH  -->  H2O+e), 
    (@reaction 0.0000000005 , C_m + H  -->  CH + e), 
    (@reaction 0.0000000000001 , C_m + H2  -->  CH2 + e), 
    (@reaction 0.0000000005 , C_m + O  -->  CO + e), 
    (@reaction 0.0000000005 , O_m + H  -->  OH + e), 
    (@reaction 0.0000000007 , O_m + H2  -->  H2O + e), 
    (@reaction 0.0000000005 , O_m + C  -->  CO + e), 
    (@reaction 0.0000000000000001 , H2+ H_p -->  H3_p ), 
    (@reaction 0.00000000000000225 , C+e  -->  C_m ), 
    (@reaction 0.00000000000000001 , C+H  -->  CH ), 
    (@reaction 0.00000000000000001 , C+H2  -->  CH2 ), 
    (@reaction 4360000000000000000 * (T/300)^(0.35) * exp(-161.3/T) , C+C  -->  C2 ), 
    (@reaction 2.1e-19 , C + O  --> CO ), 
    (@reaction 0.000000000000000446 * T^(-0.5) * exp(-4.93 / (T^(2/3)) ) , C_p + H  -->  CH_p  ), 
    (@reaction 0.0000000000000004 *(T/300)^(-0.2) , C_p + H2  -->  CH2_p ), 
    (@reaction 2.5e-18 , C_p + O  -->  CO_p  ), 
    (@reaction 0.0000000000000015 , O+e  -->  O_m ), 
    (@reaction 9.9e-19*(T/300)^(-0.38) , O+H  -->  OH ), 
    (@reaction 4.9e-20 * (T/300)^1.58 , O+O  -->  O2 ), 
    (@reaction 5.26e-18 *(T/300)^(-5.22) *exp(-90/T) , OH+H  -->  H2O ), 
    (@reaction 1.32e-32 * (T/300)^-0.38 , H+H+H  -->  H2+H), 
    (@reaction 2.8e-31 * (T)^-0.6 , H+H+H2  -->  H2+H2), 
    (@reaction 6.9e-32 * (T)^(-0.4) , H+H+He  -->  H2+He), 
    (@reaction 5.99e-33 * (T/5000)^(-1.6) , C+C+M  -->  C2+M), 
    (@reaction 6.16e-29 * (T/300)^(-3.08) , C+O+M  -->  CO+M), 
    #(@reaction 100 * C210 , C_p + O+M  -->  CO_p + M), 
    #(@reaction 100 * C210 , C+O_p + M  -->  CO_p + M), 
    (@reaction 4.33e-32 * (T/300)^(-1) , O+H+M  -->  OH+M), 
    (@reaction 2.56e-31 * (T/300)^(-2) , OH+H+M  -->  H2O+M), 
    (@reaction 9.2e-34 * (T/300)^(-1) , O+O+M  -->  O2+M), 
    (@reaction 0.00000000002 * (T/300)^(0.44) , O+CH  -->  HCO_p + e), 
    #(@reaction 0.000000000000000003 * (T)^(0.5) * (1 + 100000 * exp(-600/T)) * (1 + 0.04*(T + T)^(0.5) + 0.002*T + 0.000008 * T^2)^(-1) , H+H(s)  -->  H2), 
    (@reaction 0.00000071 , H_m    -->  H+e), 
    (@reaction 0.0000000011 , H2_p   -->  H+H_p), 
    (@reaction 0.000000000056 , H2   -->  H+H), 
    (@reaction 0.00000000000049 , H3_p   -->  H2+H_p), 
    (@reaction 0.00000000000049 , H3_p   -->  H2_p+H), 
    (@reaction 0.00000000031 , C   -->  C_p + e), 
    (@reaction 0.00000024 , C_m  -->  C+e), 
    (@reaction 0.00000000087 , CH   -->  C+H), 
    (@reaction 0.00000000077 , CH   -->  CH_p + e), 
    (@reaction 0.00000000026 , CH_p    -->  C+H_p), 
    (@reaction 0.00000000071 , CH2   -->  CH+H), 
    (@reaction 0.00000000059 , CH2   -->  CH2_p+e), 
    (@reaction 0.00000000046 , CH2_p   -->  CH_p + H), 
    (@reaction 0.000000001 , CH3_p   -->  CH2_p+H), 
    (@reaction 0.000000001 , CH3_p   -->  CH_p + H2), 
    (@reaction 0.00000000015 , C2   -->  C+C), 
    (@reaction 0.00000024 , O_m -->  O+e), 
    (@reaction 0.00000000037 , OH   -->  O+H), 
    (@reaction 0.0000000000016 , OH   -->  OH_p + e), 
    (@reaction 0.000000000001 , OH_p    -->  O+H_p), 
    (@reaction 0.0000000006 , H2O   -->  OH+H), 
    (@reaction 0.000000000032 , H2O   -->  H2O_p + e), 
    (@reaction 0.00000000005 , H2O_p    -->  H2_p+O), 
    (@reaction 0.00000000005 , H2O_p    -->  H_p + OH), 
    (@reaction 0.00000000005 , H2O_p    -->  O_p + H2), 
    (@reaction 0.00000000015 , H2O_p    -->  OH_p + H), 
    (@reaction 0.000000000025 , H3O_p    -->  H_p + H2O), 
    (@reaction 0.000000000025 , H3O_p    -->  H2_p+OH), 
    (@reaction 0.0000000000075 , H3O_p    -->  H2O_p + H), 
    (@reaction 0.000000000025 , H3O_p    -->  OH_p + H2), 
    (@reaction 0.000000000056 , O2   -->  O2_p+e), 
    (@reaction 0.0000000007 , O2   -->  O+O), 
    (@reaction 0.0000000002 , CO   -->  C+O), 
    (@reaction 1 , H   -->  H_p + e), 
    (@reaction 1.1 , He   -->  He_p + e), 
    (@reaction 0.037 , H2   -->  H_p + H+e), 
    (@reaction 0.22 , H2   -->  H+H), 
    (@reaction 0.00065 , H2   -->  H_p), 
    (@reaction 2 , H2   -->  H2_p+e), 
    (@reaction 3.8 , C   -->  C_p + e), 
    (@reaction 5.7 , O   -->  O_p + e), 
    (@reaction 6.5 , CO   -->  CO_p + e), 
    (@reaction 2800 , C   -->  C_p + e), 
    (@reaction 4000 , CH   -->  C+H), 
    (@reaction 960 , CH_p    -->  C_p + H), 
    (@reaction 2700 , CH2   -->  CH2_p+e), 
    (@reaction 2700 , CH2   -->  CH+H), 
    (@reaction 1300 , C2   -->  C+C), 
    (@reaction 2800 , OH   -->  O+H), 
    (@reaction 5300 , H2O   -->  OH+H), 
    (@reaction 4100 , O2   -->  O+O), 
    (@reaction 640 , O2   -->  O2_p+e), 
    #(@reaction 0.21*T^(1/2) *1 *(0.00001)^-1/2 , CO   -->  C+O)
    
]


    # %% Turn the network into an ODE system
@named system = ReactionSystem(reaction_equations, t)
odesys = convert(ODESystem, complete(system))
sys = convert(ODESystem, complete(system))
ssys = structural_simplify(sys)

end