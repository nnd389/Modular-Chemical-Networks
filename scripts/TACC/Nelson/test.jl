using Catalyst
using DifferentialEquations
using Random

# Define the reaction network
@variables t
@species A(t) B(t)
@parameters k1 k2

rxs = [Reaction(k1, [A], [B]),   # A → B
       Reaction(k2, [B], [A])]   # B → A
#CUDA.allowscalar(false)

# Set initial conditions and parameters
u0 = [1.0, 0.0]
p = [1.0, 0.5]
tspan = (0.0, 10.0)

@named system = ReactionSystem(rxs, t)
odesys = Catalyst.complete(system) # Both Catalyst and CUDA uhave their own "complete" so must specify
prob = ODEProblem(odesys, u0, tspan, p)
sol = solve(prob)
print("\nRegular Solved!")


u0_gpu = Float32[1.0, 0.0]
p_gpu = Float32[1.0, 0.5]
tspan_gpu = (0.0f0, 10.0f0)

Random.seed!(1)
function prob_func_distributed(prob, i, repeat)
        remake(prob, u0 = rand(2) .* u0_gpu)
end

prob_ens = ODEProblem(odesys, u0_gpu, tspan_gpu, p_gpu)
monteprob_simple = EnsembleProblem(prob_ens, prob_func = prob_func_distributed, safetycopy = false);
sol_ensemble = solve(monteprob_simple, EnsembleThreads(), trajectories = 10)
print("\nEnsemble Solved!")
using CUDA
using DiffEqGPU
sol_gpu = solve(monteprob_simple, Rodas4(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = 10)
print("\nCUDA Solved!\n")