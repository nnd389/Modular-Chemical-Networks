

### Let's start testing out the documentation from ###
### Setting Up Multi-GPU Parallel Parameter Sweeps ###

print("\n\n GPU Batch Problem")
num_runs = 1000

# Setup processes with different CUDA devices
using DiffEqGPU
using OrdinaryDiffEq
using Test
using Distributed
using Test
using Random
using SparseArrays


print("\nGPU Batch 1: Finished the packages")



@everywhere begin
        function Nelson_gpu(du,u,p,t)
                T, Av, Go, n_H, shield = p
            # 1: H2
            #@inbounds begin
            du[1] = -1.2f-17 * u[1] + 
                    n_H * (1.9f-6 * u[2] * u[3]) / (T^0.54) - 
                    n_H * 4f-16 * u[1] * u[12] - 
                    n_H * 7f-15 * u[1] * u[5] + 
                    n_H * 1.7f-9 * u[10] * u[2] + 
                    n_H * 2f-9 * u[2] * u[6] + 
                    n_H * 2f-9 * u[2] * u[14] + 
                    n_H * 8f-10 * u[2] * u[8] 
            
            # 2: H3+
            du[2] = 1.2f-17 * u[1] + 
                    n_H * (-1.9f-6 * u[3] * u[2]) / (T^0.54) - 
                    n_H * 1.7f-9 * u[10] * u[2] - 
                    n_H * 2f-9 * u[2] * u[6] - 
                    n_H * 2f-9 * u[2] * u[14] - 
                    n_H * 8f-10 * u[2] * u[8]
            
            # 3: e
            du[3] = n_H * (-1.4f-10 * u[3] * u[12]) / (T^0.61) - 
                    n_H * (3.8f-10 * u[13] * u[3]) / (T^0.65) - 
                    n_H * (3.3f-5 * u[11] * u[3]) / T + 
                    1.2f-17 * u[1] - 
                    n_H * (1.9f-6 * u[3] * u[2]) / (T^0.54) + 
                    6.8f-18 * u[4] - 
                    n_H * (9f-11 * u[3] * u[5]) / (T^0.64) + 
                    3f-10 * Go * exp(-3 * Av) * u[6] +
                    n_H * 2f-9 * u[2] * u[13] # added this extra term from a CR ionization reaction
                    + 2.0f-10 * Go * exp(-1.9 * Av) * u[14] # this term was added as part of the skipped photoreaction
            
            
            # 4: He
            du[4] = n_H * (9f-11 * u[3] * u[5]) / (T^0.64) - 
                    6.8f-18 * u[4] + 
                    n_H * 7f-15 * u[1] * u[5] + 
                    n_H * 1.6f-9 * u[10] * u[5]
            
            # 5: He+
            du[5] = 6.8f-18 * u[4] - 
                    n_H * (9f-11 * u[3] * u[5]) / (T^0.64) - 
                    n_H * 7f-15 * u[1] * u[5] - 
                    n_H * 1.6f-9 * u[10] * u[5]
            
            # 6: C
            du[6] = n_H * (1.4f-10 * u[3] * u[12]) / (T^0.61) - 
                    n_H * 2f-9 * u[2] * u[6] - 
                    n_H * 5.8f-12 * (T^0.5) * u[9] * u[6] + 
                    1f-9 * Go * exp(-1.5 * Av) * u[7] - 
                    3f-10 * Go * exp(-3 * Av) * u[6] + 
                    1f-10 * Go * exp(-3 * Av) * u[10] * shield
            
            # 7: CHx
            du[7] = n_H * (-2f-10) * u[7] * u[8] + 
                    n_H * 4f-16 * u[1] * u[12] + 
                    n_H * 2f-9 * u[2] * u[6] - 
                    1f-9 * Go * u[7] * exp(-1.5 * Av)
            
            # 8: O
            du[8] = n_H * (-2f-10) * u[7] * u[8] + 
                    n_H * 1.6f-9 * u[10] * u[5] - 
                    n_H * 8f-10 * u[2] * u[8] + 
                    5f-10 * Go * exp(-1.7 * Av) * u[9] + 
                    1f-10 * Go * exp(-3 * Av) * u[10] * shield
            
            # 9: OHx
            du[9] = n_H * (-1f-9) * u[9] * u[12] + 
                    n_H * 8f-10 * u[2] * u[8] - 
                    n_H * 5.8f-12 * (T^0.5) * u[9] * u[6] - 
                    5f-10 * Go * exp(-1.7 * Av) * u[9]
            
            # 10: CO
            du[10] = n_H * (3.3f-5 * u[11] * u[3]) / T + 
                    n_H * 2f-10 * u[7] * u[8] - 
                    n_H * 1.7f-9 * u[10] * u[2] - 
                    n_H * 1.6f-9 * u[10] * u[5] + 
                    n_H * 5.8f-12 * (T^0.5) * u[9] * u[6] - 
                    1f-10 * Go * exp(-3 * Av) * u[10] + 
                    1.5f-10 * Go * exp(-2.5 * Av) * u[11] * shield
            
            # 11: HCO+
            du[11] = n_H * (-3.3f-5 * u[11] * u[3]) / T + 
                    n_H * 1f-9 * u[9] * u[12] + 
                    n_H * 1.7f-9 * u[10] * u[2] - 
                    1.5f-10 * Go * exp(-2.5 * Av) * u[11]
            
            # 12: C+
            du[12] = n_H * (-1.4f-10 * u[3] * u[12]) / (T^0.61) - 
                    n_H * 4f-16 * u[1] * u[12] - 
                    n_H * 1f-9 * u[9] * u[12] + 
                    n_H * 1.6f-9 * u[10] * u[5] + 
                    3f-10 * Go * exp(-3 * Av) * u[6]
            
            # 13: M+
            du[13] = n_H * (-3.8f-10 * u[13] * u[3]) / (T^0.65) + 
                    n_H * 2f-9 * u[2] * u[14] 
                    + 2.0f-10 * Go * exp(-1.9 * Av) * u[14] # this term was added as part of the skipped photoreaction
            
            # 14: M
            du[14] = n_H * (3.8f-10 * u[13] * u[3]) / (T^0.65) - 
                    n_H * 2f-9 * u[2] * u[14] 
                    - 2.0f-10 * Go * exp(-1.9 * Av) * u[14] # this term was added as part of the skipped photoreaction
            #end
        end
 
        #J = spzeros(Float32, 14, 14) # This is for sparse matrices, I hope to use it one day

        function Nelson_Jac(J, u, p, t)
                T, Av, Go, n_H, shield = p
            
                fill!(J, 0.0f0)
                #fill!(J.nzval, 0.0f0)
                
            
                # du[1]
                J[1,1] = -1.2f-17 - n_H * 4f-16 * u[12] - n_H * 7f-15 * u[5]
                J[1,2] = n_H * (1.9f-6 * u[3]) / (T^0.54) + n_H * 1.7f-9 * u[10] + n_H * 2f-9 * u[6] + n_H * 2f-9 * u[14] + n_H * 8f-10 * u[8]
                J[1,3] = n_H * (1.9f-6 * u[2]) / (T^0.54)
                J[1,5] = -n_H * 7f-15 * u[1]
                J[1,6] = n_H * 2f-9 * u[2]
                J[1,8] = n_H * 8f-10 * u[2]
                J[1,10] = n_H * 1.7f-9 * u[2]
                J[1,12] = -n_H * 4f-16 * u[1]
                J[1,14] = n_H * 2f-9 * u[2]
            
                # du[2]
                J[2,1] = 1.2f-17
                J[2,2] = -n_H * (1.9f-6 * u[3]) / (T^0.54) - n_H * 1.7f-9 * u[10] - n_H * 2f-9 * u[6] - n_H * 2f-9 * u[14] - n_H * 8f-10 * u[8]
                J[2,3] = -n_H * (1.9f-6 * u[2]) / (T^0.54)
                J[2,6] = -n_H * 2f-9 * u[2]
                J[2,8] = -n_H * 8f-10 * u[2]
                J[2,10] = -n_H * 1.7f-9 * u[2]
                J[2,14] = -n_H * 2f-9 * u[2]
            
                # du[3]
                J[3,2] = -n_H * (1.9f-6 * u[3]) / (T^0.54) + n_H * 2f-9 * u[13]
                J[3,3] = -n_H * (1.4f-10 * u[12]) / (T^0.61) - n_H * (3.8f-10 * u[13]) / (T^0.65) - n_H * (3.3f-5 * u[11]) / T - n_H * (1.9f-6 * u[2]) / (T^0.54) - n_H * (9f-11 * u[5]) / (T^0.64)
                J[3,4] = 6.8f-18
                J[3,5] = -n_H * (9f-11 * u[3]) / (T^0.64)
                J[3,6] = 3f-10 * Go * exp(-3 * Av)
                J[3,8] = 0.0
                J[3,11] = -n_H * (3.3f-5 * u[3]) / T
                J[3,12] = -n_H * (1.4f-10 * u[3]) / (T^0.61)
                J[3,13] = -n_H * (3.8f-10 * u[3]) / (T^0.65)
                J[3,14] = 2.0f-10 * Go * exp(-1.9 * Av)
            
                # du[4]
                J[4,1] = n_H * 7f-15 * u[5]
                J[4,3] = n_H * (9f-11 * u[5]) / (T^0.64)
                J[4,4] = -6.8f-18
                J[4,5] = n_H * (9f-11 * u[3]) / (T^0.64) + n_H * 1.6f-9 * u[10]
                J[4,10] = n_H * 1.6f-9 * u[5]
            
                # du[5]
                J[5,1] = -n_H * 7f-15 * u[5]
                J[5,3] = -n_H * (9f-11 * u[5]) / (T^0.64)
                J[5,4] = 6.8f-18
                J[5,5] = -n_H * (9f-11 * u[3]) / (T^0.64) - n_H * 7f-15 * u[1] - n_H * 1.6f-9 * u[10]
                J[5,10] = -n_H * 1.6f-9 * u[5]
            
                # du[6]
                J[6,2] = -n_H * 2f-9 * u[2]
                J[6,3] = n_H * (1.4f-10 * u[12]) / (T^0.61)
                J[6,6] = -n_H * 2f-9 * u[2] - n_H * 5.8f-12 * (T^0.5) * u[9] - 3f-10 * Go * exp(-3 * Av)
                J[6,7] = 0.0
                J[6,9] = -n_H * 5.8f-12 * (T^0.5) * u[6]
                J[6,10] = 1f-10 * Go * exp(-3 * Av) * shield
                J[6,12] = n_H * (1.4f-10 * u[3]) / (T^0.61)
            
                # du[7]
                J[7,1] = n_H * 4f-16 * u[12]
                J[7,2] = n_H * 2f-9 * u[6]
                J[7,6] = n_H * 2f-9 * u[2]
                J[7,7] = -n_H * 2f-10 * u[8] - 1f-9 * Go * exp(-1.5 * Av)
                J[7,8] = -n_H * 2f-10 * u[7]
                J[7,12] = n_H * 4f-16 * u[1]
            
                # du[8]
                J[8,2] = -n_H * 8f-10 * u[2]
                J[8,5] = n_H * 1.6f-9 * u[10]
                J[8,7] = -n_H * 2f-10 * u[8]
                J[8,8] = -n_H * 2f-10 * u[7] - n_H * 8f-10 * u[2]
                J[8,9] = 5f-10 * Go * exp(-1.7 * Av)
                J[8,10] = 1f-10 * Go * exp(-3 * Av) * shield
                J[8,10] += n_H * 1.6f-9 * u[5]
            
                # du[9]
                J[9,2] = n_H * 8f-10 * u[8]
                J[9,6] = -n_H * 5.8f-12 * (T^0.5) * u[9]
                J[9,8] = n_H * 8f-10 * u[2]
                J[9,9] = -n_H * 1f-9 * u[12] - n_H * 5.8f-12 * (T^0.5) * u[6] - 5f-10 * Go * exp(-1.7 * Av)
                J[9,12] = -n_H * 1f-9 * u[9]
            
                # du[10]
                J[10,2] = -n_H * 1.7f-9 * u[10]
                J[10,3] = n_H * (3.3f-5 * u[11]) / T
                J[10,5] = -n_H * 1.6f-9 * u[10]
                J[10,6] = n_H * 5.8f-12 * (T^0.5) * u[9]
                J[10,7] = n_H * 2f-10 * u[8]
                J[10,8] = n_H * 2f-10 * u[7]
                J[10,9] = n_H * 5.8f-12 * (T^0.5) * u[6]
                J[10,10] = -n_H * 1.7f-9 * u[2] - n_H * 1.6f-9 * u[5] - 1f-10 * Go * exp(-3 * Av)
                J[10,11] = n_H * (3.3f-5 * u[3]) / T + 1.5f-10 * Go * exp(-2.5 * Av) * shield
            
                # du[11]
                J[11,2] = n_H * 1.7f-9 * u[10]
                J[11,3] = -n_H * (3.3f-5 * u[11]) / T
                J[11,9] = n_H * 1f-9 * u[12]
                J[11,10] = n_H * 1.7f-9 * u[2]
                J[11,11] = -n_H * (3.3f-5 * u[3]) / T - 1.5f-10 * Go * exp(-2.5 * Av)
                J[11,12] = n_H * 1f-9 * u[9]
            
                # du[12]
                J[12,1] = -n_H * 4f-16 * u[12]
                J[12,3] = -n_H * (1.4f-10 * u[12]) / (T^0.61)
                J[12,5] = n_H * 1.6f-9 * u[10]
                J[12,6] = 3f-10 * Go * exp(-3 * Av)
                J[12,9] = -n_H * 1f-9 * u[12]
                J[12,10] = n_H * 1.6f-9 * u[5]
                J[12,12] = -n_H * (1.4f-10 * u[3]) / (T^0.61) - n_H * 4f-16 * u[1] - n_H * 1f-9 * u[9]
            
                # du[13]
                J[13,2] = n_H * 2f-9 * u[14]
                J[13,3] = -n_H * (3.8f-10 * u[13]) / (T^0.65)
                J[13,13] = -n_H * (3.8f-10 * u[3]) / (T^0.65)
                J[13,14] = n_H * 2f-9 * u[2] + 2.0f-10 * Go * exp(-1.9 * Av)
            
                # du[14]
                J[14,2] = -n_H * 2f-9 * u[14]
                J[14,3] = n_H * (3.8f-10 * u[13]) / (T^0.65)
                J[14,13] = n_H * (3.8f-10 * u[3]) / (T^0.65)
                J[14,14] = -n_H * 2f-9 * u[2] - 2.0f-10 * Go * exp(-1.9 * Av)
            end










        u0_gpu = Float32[0.5;    # 1:  H2   yep?
        9.059f-9; # 2:  H3+  yep
        2.0f-4;   # 3:  e    yep
        0.1;      # 4:  He  SEE lines 535 NL99
        7.866f-7; # 5:  He+  yep? should be 2.622f-5
        0.0;      # 6:  C    yep
        0.0;      # 7:  CHx  yep
        0.0004;   # 8:  O    yep
        0.0;      # 9:  OHx  yep
        0.0;      # 10: CO   yep
        0.0;      # 11: HCO+ yep
        0.0002;   # 12: C+   yep
        2.0f-7;   # 13: M+   yep
        2.0f-7]   # 14: M    yep

        tspan_gpu = (0.0f0, 946080000000.0f0) # ~30 thousand yrs
        params_gpu = Float32[10.0f0, 2.0f0, 1.7f0, 611.0f0, 1.0f0]

        Random.seed!(1)
        function prob_func_distributed(prob_gpu, i, repeat)
            remake(prob_gpu, u0 = rand(14) .* u0_gpu)
        end


            print("\nGPU Batch 3: Finished @everywhere macro")
            print("\nTypes!!")
            print("\nType of u0_gpu: ", typeof(u0_gpu))
            print("\nType of tspan_gpu: ", typeof(tspan_gpu))
            print("\nType of params_gpu: ", typeof(params_gpu))



end

f_gpu = ODEFunction(Nelson_gpu, jac = Nelson_Jac)
prob_gpu = ODEProblem(f_gpu, u0_gpu, tspan_gpu, params_gpu)


print("\nGPU Batch: 4: Finished spawnat")
#prob_gpu = ODEProblem(Nelson_gpu, u0_gpu, tspan_gpu, params_gpu)
print("\nGPU Batch 5: Finished ODEProblem")
monteprob_gpu = EnsembleProblem(prob_gpu, prob_func = prob_func_distributed, safetycopy = false) # try taking out safeteycopy=false
print("\nGPU Batch 6: Finished Ensemble Problem")


# NOTE: documentation examples usually randomize the parameters instead of initial conditions
print("\n We are about to solve the gpu batch problem with ", num_runs, " runs and batch size of ", num_runs/2)



print("\n\nGPU batch: GPU Timing to solve ", num_runs, " random Nelson systems on GPUs with Rodas5P():")
@time solve(monteprob_gpu, Rodas5P(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)
@time solve(monteprob_gpu, Rodas5P(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)
@time solve(monteprob_gpu, Rodas5P(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)
@time solve(monteprob_gpu, Rodas5P(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)


print("\n\nGPU batch: GPU Timing to solve ", num_runs, " random Nelson systems on GPUs with Rodas4():")
@time solve(monteprob_gpu, Rodas4(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)
@time solve(monteprob_gpu, Rodas4(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)
@time solve(monteprob_gpu, Rodas4(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)
@time solve(monteprob_gpu, Rodas4(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)


#print("\n\nGPU batch: GPU Timing to solve ", num_runs, " random Nelson systems on GPUs with Tsit5():")
#print("\n first solve")
#@time solve(monteprob_gpu, Tsit5(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)
#@time solve(monteprob_gpu, Tsit5(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)


#print("\n\nGPU batch: GPU Timing to solve ", num_runs, " random Nelson systems on GPUs with Vern9():")
#@time solve(monteprob_gpu, Vern9(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)
#@time solve(monteprob_gpu, Vern9(), EnsembleThreads(), trajectories = num_runs, batch_size = num_runs, reltol=1.49012f-8, abstol=1.49012f-8, saveat=1e10)


