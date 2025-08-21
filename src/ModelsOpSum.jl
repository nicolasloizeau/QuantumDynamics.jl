module ModelsOpSum

using ITensors, ITensorMPS

"""
    ising(N, hx, hz)

Hamiltonian for the 1D Ising model with periodic boundary conditions.
``H= \\sum_{j=1}^N ( Z_j Z_{j+1} + h_x X_j + h_z Z_j)``
"""
function ising(N, hx, hz)
    H = OpSum()
    for j in 1:N
        H += 4 , "Sz", j, "Sz", mod1(j + 1, N)
        H += 2*hx, "Sx", j
        H += 2*hz, "Sz", j
    end
    return H
end


function XYZ(N, alpha)
    H = OpSum()
    for j in 1:N
        H += 2, alpha, j
    end
    return H
end

X(N) = XYZ(N, "Sx")
Y(N) = XYZ(N, "Sy")
Z(N) = XYZ(N, "Sz")

end
