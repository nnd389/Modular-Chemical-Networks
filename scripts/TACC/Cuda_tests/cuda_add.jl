using CUDA

# Function for vector addition on the GPU
function vector_add_gpu(A, B)
    # Copy input vectors to the GPU
    d_A = CUDA.fill(0.0f0, length(A))  # Create a device array with the same size as A
    d_B = CUDA.fill(0.0f0, length(B))  # Create a device array with the same size as B
    CUDA.copyto!(d_A, A)               # Copy data from host to device
    CUDA.copyto!(d_B, B)               # Copy data from host to device

    # Allocate an output array on the GPU
    d_C = CUDA.fill(0.0f0, length(A))

    # Launch the kernel (this is a simple vector add kernel)
    @cuda threads=256 vector_add_kernel(d_A, d_B, d_C)

    # Copy the result back to the host
    C = Array(d_C)  # Convert the device array back to a host array
    return C
end

# Define a simple kernel for vector addition (GPU version)
function vector_add_kernel(A, B, C)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x  # Global index
    if i <= length(A)
        C[i] = A[i] + B[i]
    end
end

# Example Usage
N = 1_000_000  # Size of the vectors
A = rand(Float32, N)  # Random vector A
B = rand(Float32, N)  # Random vector B

# Perform the vector addition on the GPU
C = vector_add_gpu(A, B)

# Check the result (optional, but for validation)
println("First 10 elements of the result: ", C[1:10])

