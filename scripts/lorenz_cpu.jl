#using DiffeqGPU, CUDA
using OrdinaryDiffEq

function lorenz!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3])
    du[3] = u[1] * u[2] - p[3] * u[3]
end


u0 = Float32[1.0; 0.0; 0.0]
tspan = (0.0f0, 100.0f0)
p = [10.0f0, 28.0f0, 8/3.0f0]
prob = ODEProblem(lorenz!, u0, tspan, p)
print("Time to solve Lorenz ONCE: ")
@time solve(prob, Tsit5(), saveat = 1.0)
@time solve(prob, Tsit5(), saveat = 1.0)
@time solve(prob, Tsit5(), saveat = 1.0)
@time solve(prob, Tsit5(), saveat = 1.0)

traj_num = 10000
prob_func = (prob, i, repeat) -> remake(prob, p=rand(Float32,3) .* p);
monteprob = EnsembleProblem(prob, prob_func = prob_func, safetycopy = false);
sol = solve(monteprob, Tsit5(), EnsembleThreads(), trajectories = traj_num, saveat = 1.0f0);
print("Ensemble Timing to solve ", traj_num, " random Lorenz systems:")
@time solve(monteprob, Tsit5(), EnsembleThreads(), trajectories = traj_num, saveat = 1.0f0);
@time solve(monteprob, Tsit5(), EnsembleThreads(), trajectories = traj_num, saveat = 1.0f0);
@time solve(monteprob, Tsit5(), EnsembleThreads(), trajectories = traj_num, saveat = 1.0f0);
@time solve(monteprob, Tsit5(), EnsembleThreads(), trajectories = traj_num, saveat = 1.0f0);


print("For Loop Timing to solve ", traj_num, " random Lorenz systems")
@time begin
    for i in 1:traj_num
        p_rand_lorenz = rand(Float32,3) .* p
        prob_rand_lorenz = ODEProblem(lorenz!, u0, tspan, p)
        sol_rand_lorenz = solve(prob_rand_lorenz, Tsit5())
    end
end
@time begin
    for i in 1:traj_num
        p_rand_lorenz = rand(Float32,3) .* p
        prob_rand_lorenz = ODEProblem(lorenz!, u0, tspan, p)
        sol_rand_lorenz = solve(prob_rand_lorenz, Tsit5())
    end
end
@time begin
    for i in 1:traj_num
        p_rand_lorenz = rand(Float32,3) .* p
        prob_rand_lorenz = ODEProblem(lorenz!, u0, tspan, p)
        sol_rand_lorenz = solve(prob_rand_lorenz, Tsit5())
    end
end
@time begin
    for i in 1:traj_num
        p_rand_lorenz = rand(Float32,3) .* p
        prob_rand_lorenz = ODEProblem(lorenz!, u0, tspan, p)
        sol_rand_lorenz = solve(prob_rand_lorenz, Tsit5())
    end
end



#plot(sol)