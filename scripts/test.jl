

function random_vector()
    vec = ones(14) 
    #u0 = [H2, H3+, e,  He, He+, C, CHx, O,  OHx, CO, HCO+, C+, M+, M]
    # I am following T. Grassi's CODE as close as possible, 
    # But it doesn't make sense to me with what the paper says
    ngas = 1e4
    rmax = -4
    rmin = -6

    vec[1] = ngas         # H2
    vec[2] = 1e-20 * ngas # H3+
    vec[3] = 4*(1e-20 * ngas) + ngas * 10^(rand() * (rmax-rmin) + rmin) # e = H3+ & He+ & HCO+ & C+ & M+ 
    vec[4] = 1e-20 * ngas # He
    vec[5] = 1e-20 * ngas # He+
    vec[6] = ngas * 10^(rand() * (rmax-rmin) + rmin) # C 
    vec[7] = 1e-20 * ngas # CHx
    vec[8] = 4*(ngas * 10^(rand() * (rmax-rmin) + rmin)) # O
    vec[9] = 1e-20 * ngas # OHx
    vec[10] = 1e-20 * ngas # CO
    vec[11] = 1e-20 * ngas # HCO+
    vec[12] = ngas * 10^(rand() * (rmax-rmin) + rmin)# C+
    vec[13] = 1e-20 * ngas # M+
    vec[14] = 1e-20 * ngas # M
    return vec
end

# Example usage:
vec = random_vector()
println(vec)  # Prints the generated vector