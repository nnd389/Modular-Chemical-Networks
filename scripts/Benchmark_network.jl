	# %% Some basic astrochemistry constants:
kboltzmann = 1.38064852e-16  # erg / K
pmass = 1.6726219e-24  # g
dust2gas = 1e-2 # ratio
mu = 2.34
gamma_ad = 1.4
gnot = 1e1

    # %% Set The timespan, paramters, and initial conditions
	seconds_per_year = 3600 * 24 * 365
	number_density = 1e5
	minimum_fractional_density = 1e-30 * number_density
	
	tspan = (0.0, 1e6*seconds_per_year)
	params = Dict(
		:T => 100,
		:dust2gas => 0.01,
		:radiation_field => 0.1,
		:cosmic_ionisation_rate => 1e-17)
	
	u0 = [
		number_density,                  # 1: H2
		number_density * 2e-4,           # 2: O
		number_density * 1e-4,           # 3: C
		minimum_fractional_density,      # 4: O⁺
		minimum_fractional_density,      # 5: OH⁺
		minimum_fractional_density,      # 6: H
		minimum_fractional_density,      # 7: H2O⁺
		minimum_fractional_density,      # 8: H3O⁺
		minimum_fractional_density,      # 9: E
		minimum_fractional_density,      # 10: H2O
		minimum_fractional_density,      # 11: OH
		minimum_fractional_density,      # 12: C⁺
		minimum_fractional_density,      # 13: CO
		minimum_fractional_density,      # 14: CO⁺
		minimum_fractional_density,      # 15: H⁺
		minimum_fractional_density,      # 16: HCO⁺
	]
	

# CONTINUE HERE
# Try this: https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/constraint_equations/#Coupling-ODE-constraints-via-directly-building-a-ReactionSystem

ka_reaction(Tgas, α=1.0, β=1.0, γ=0.0) = α*(Tgas/300)^β*exp(−γ / Tgas)
@variables t 
@species H2(t) O(t) C(t) O⁺(t) OH⁺(t) H(t) H2O⁺(t) H3O⁺(t) E(t) H2O(t) OH(t) C⁺(t) CO(t) CO⁺(t) H⁺(t) HCO⁺(t)
@parameters T cosmic_ionisation_rate radiation_field dust2gas

D = Differential(t)
reaction_equations = [
	(@reaction 1.6e-9, O⁺ + H2 --> OH⁺ + H),
	(@reaction 1e-9, OH⁺ + H2 --> H2O⁺ + H),
	(@reaction 6.1e-10, H2O⁺ + H2 --> H3O⁺ + H),
	(@reaction ka_reaction(T, 1.1e-7, -1/2), H3O⁺ + E --> H2O + H),
	(@reaction ka_reaction(T, 8.6e-8, -1/2), H2O⁺ + E --> OH + H),
	(@reaction ka_reaction(T, 3.9e-8, -1/2), H2O⁺ + E --> O + H2),
	(@reaction ka_reaction(T, 6.3e-9, -0.48), OH⁺ + E --> O + H),
	(@reaction ka_reaction(T, 3.4e-12, -0.63), O⁺ + E --> O),
	(@reaction 2.8 * cosmic_ionisation_rate, O --> O⁺ + E),
	(@reaction 2.62 * cosmic_ionisation_rate, C --> C⁺ + E),
	(@reaction 5.0 * cosmic_ionisation_rate, CO --> C + O),
	(@reaction ka_reaction(T, 4.4e-12, -0.61), C⁺ + E --> C),
	(@reaction ka_reaction(T, 1.15e-10, -0.339), C⁺ + OH --> CO + H),
	(@reaction 9.15e-10 * (0.62 + 0.4767 * 5.5 * sqrt(300 / T)), C⁺ + OH --> CO⁺ + H),
	(@reaction 4e-10, CO⁺ + H --> CO + H⁺),
	(@reaction 7.28e-10, CO⁺ + H2 --> HCO⁺ + H),
	(@reaction ka_reaction(T, 2.8e-7, -0.69), HCO⁺ + E --> CO + H),
	(@reaction ka_reaction(T, 3.5e-12, -0.7), H⁺ + E --> H),
	(@reaction 2.121e-17 * dust2gas / 1e-2, H + H --> H2),
    (@reaction 1e-1 * cosmic_ionisation_rate, H2 --> H + H),
	(@reaction 3.39e-10 * radiation_field, C --> C⁺ + E),
	(@reaction 2.43e-10 * radiation_field, CO --> C + O),
	(@reaction 7.72e-10 * radiation_field, H2O --> OH + H),
]

    # %% Turn the network into an ODE system
@named system = ReactionSystem(reaction_equations, t)
odesys = convert(ODESystem, complete(system))
sys = convert(ODESystem, complete(system))
ssys = structural_simplify(sys) 

