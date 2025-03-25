using OrdinaryDiffEq, CUDA, LinearAlgebra

# A: 1000 square system
u0_A = cu(rand(1000));
A = cu(randn(1000, 1000));
f_A(du, u, p, t) = mul!(du, A, u);
prob_A = ODEProblem(f_A, u0_A, (0.0f0, 1.0f0)); # Float32 is better on GPUs!
print("Time to solve 1000 square system: ")
@time solve(prob_A, Tsit5());

# B: 3000 square system
u0_B = cu(rand(3000));
B = cu(randn(3000, 3000));
f_B(du, u, p, t) = mul!(du, B, u);
prob_B = ODEProblem(f_B, u0_B, (0.0f0, 1.0f0)); # Float32 is better on GPUs!
print("Time to solve 3000 square system: ")
@time solve(prob_B, Tsit5());

# C: 5000 square system
u0_C = cu(rand(5000));
C = cu(randn(5000, 5000));
f_C(du, u, p, t) = mul!(du, C, u);
prob_C = ODEProblem(f_C, u0_C, (0.0f0, 1.0f0)); # Float32 is better on GPUs!
print("Time to solve 5000 square system: ")
@time solve(prob_C, Tsit5());
