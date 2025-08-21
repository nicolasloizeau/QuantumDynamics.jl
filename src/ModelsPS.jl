
module ModelsPS

using PauliStrings

"""
    ising_ts(N, hx, hz)

Hamiltonian for the 1D Ising model with periodic boundary conditions.
``H= \\sum_{j=1}^N ( Z_j Z_{j+1} + h_x X_j + h_z Z_j)``
Returns a `OperatorTS1D` object.
"""
function ising_ts(N, hx, hz)
    H = Operator(N)
    H += "Z", 1, "Z", 2
    H += hx, "X", 1
    H += hz, "Z", 1
    return OperatorTS1D(H, full=false)
end

ising(N, hx, hz) = Operator(ising_ts(N, hx, hz))

function XYZ_ts(N, alpha)
    H = Operator(N)
    H += alpha, 1
    return OperatorTS1D(H, full=false)
end

X_ts(N) = XYZ_ts(N, "X")
Y_ts(N) = XYZ_ts(N, "Y")
Z_ts(N) = XYZ_ts(N, "Z")

"""
    X(N), Y(N), Z(N)
Total X, Y, Z operators for a system of size `N`.
"""
X(N) = Operator(X_ts(N))
Y(N) = Operator(Y_ts(N))
Z(N) = Operator(Z_ts(N))



end
