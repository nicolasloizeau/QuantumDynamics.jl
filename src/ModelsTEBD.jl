
module ModelsTEBD

using ITensors, ITensorMPS

"""
    ising(N, hx, hz, sites)

Hamiltonian for the 1D Ising model with periodic boundary conditions.
``H= \\sum_{j=1}^N ( Z_j Z_{j+1} + h_x X_j + h_z Z_j)``
"""
function ising(N, hx, hz, sites)
    terms = ITensor[]
    for j in 1:N
        s1 = sites[j]
        s2 = sites[mod1(j + 1, N)]
        push!(terms, 4 * op("Sz", s1) * op("Sz", s2))
        push!(terms, 2 * hx * op("Sx", s1))
        push!(terms, 2 * hz * op("Sz", s1))
    end
    return terms
end


end
