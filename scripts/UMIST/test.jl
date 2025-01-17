# This will find all the perumations of a 0 and 1 array of length three
#=
# Loop through all combinations of replacing elements with 1
for i in 0:7  # This gives us 8 combinations (from 0b000 to 0b111)
    # Create a copy of the original array
    u = copy(u0_run)
    
    # Use a bitwise operation to set the corresponding positions to 1
    for j in 1:3
        if (i >> (j-1)) & 1 == 1
            u[j] = 1
        end
    end
    
    # Print the resulting array
    println(u)
end
=#

@time begin
for i = 1:100000
    random_array = rand(3)
end

end

