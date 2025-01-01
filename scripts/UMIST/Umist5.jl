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
u0 = [
      0.0,   # 1: C⁻
      0.0,   # 2: C
      0.0,   # 3: C2
      0.0,   # 4: e
      0.0,   # 5: CH2
      0.0,   # 6: C2H2
      0.0,   # 7: CH
      0.0,   # 8: C2H
      0.0,   # 9: CO2
      0.0,   # 10: CO
      0.0,   # 11: H2O
      0.0]   # 12: H2CO


### Parameters ###
T = 10
params = Dict(
         :T => T,
         :Av => 2,
         :n_H => 611,
         :cr_ion_rate => 6e-18,
         :omega => .5,
         :k0 => 5e-10 * (T/300)^(0) * exp(-0/T),  # Reaction rate number 0 with type AD
         :k1 => 5e-10 * (T/300)^(0) * exp(-0/T),  # Reaction rate number 1 with type AD
         :k2 => 5e-10 * (T/300)^(0) * exp(-0/T),  # Reaction rate number 2 with type AD
         :k3 => 4.7e-11 * (T/300)^(0) * exp(-0/T),  # Reaction rate number 3 with type AD
         :k4 => 5e-10 * (T/300)^(0) * exp(-0/T) # Reaction rate number 4 with type AD
         )


### Tempurature range flags ###
T_low_bounds = [10 10 10 10 10 ]
T_upp_bounds = [41000 41000 41000 41000 41000 ]
for i = 1:length(T_low_bounds)
    if T < T_low_bounds[i]
        print("Tempurature below the tempurature range [", T_low_bounds[i], ", ", T_upp_bounds[i], "] for reaction number ", i ) 
    elseif T > T_upp_bounds[i]
        print("Tempurature higher than the tempurature range [", T_low_bounds[i], ", ", T_upp_bounds[i], "] for reaction number ", i ) 
    end
end


### Network Admin things ###
@variables t 
@species C⁻(t) C(t) C2(t) e(t) CH2(t) C2H2(t) CH(t) C2H(t) CO2(t) CO(t) H2O(t) H2CO(t) 
@parameters T Av n_H cr_ion_rate omega k0 k1 k2 k3 k4 


### Reaction Equations ###
reaction_equations = [
                      (@reaction n_H * k0, C⁻ + C --> C2 + e), 
                      (@reaction n_H * k1, C⁻ + CH2 --> C2H2 + e), 
                      (@reaction n_H * k2, C⁻ + CH --> C2H + e), 
                      (@reaction n_H * k3, C⁻ + CO2 --> CO + CO + e), 
                      (@reaction n_H * k4, C⁻ + H2O --> H2CO + e)] 


### Turn the Network into a system of ODEs ###
@named system = ReactionSystem(reaction_equations, t)
odesys = convert(ODESystem, complete(system))
sys = convert(ODESystem, complete(system))
ssys = structural_simplify(sys)
prob = ODEProblem(ssys, u0, tspan, params)
sol = solve(prob, Rodas4())


### Plotting ###
# Species number 0: C⁻
plot(sol, idxs = (0,1), lw = 3, lc = "blue", title = "Abundance of C⁻")

# Species number 1: C
plot!(sol, idxs = (0,2), lw = 3, lc = "blue", title = "Abundance of C")

#= 

# Species number 2: C2
plot(sol, idxs = (0,2), lw = 3, lc = "blue", title = "Abundance of C2")

# Species number 3: e
plot(sol, idxs = (0,3), lw = 3, lc = "blue", title = "Abundance of e")

# Species number 4: CH2
plot(sol, idxs = (0,4), lw = 3, lc = "blue", title = "Abundance of CH2")

# Species number 5: C2H2
plot(sol, idxs = (0,5), lw = 3, lc = "blue", title = "Abundance of C2H2")

# Species number 6: CH
plot(sol, idxs = (0,6), lw = 3, lc = "blue", title = "Abundance of CH")

# Species number 7: C2H
plot(sol, idxs = (0,7), lw = 3, lc = "blue", title = "Abundance of C2H")

# Species number 8: CO2
plot(sol, idxs = (0,8), lw = 3, lc = "blue", title = "Abundance of CO2")

# Species number 9: CO
plot(sol, idxs = (0,9), lw = 3, lc = "blue", title = "Abundance of CO")

# Species number 10: H2O
plot(sol, idxs = (0,10), lw = 3, lc = "blue", title = "Abundance of H2O")

# Species number 11: H2CO
plot(sol, idxs = (0,11), lw = 3, lc = "blue", title = "Abundance of H2CO")
=#