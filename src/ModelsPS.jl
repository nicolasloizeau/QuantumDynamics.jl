
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

function tot_ts(N, alpha)
    H = Operator(N)
    H += alpha, 1
    return OperatorTS1D(H, full=false)
end

X_ts(N) = tot_ts(N, "X")
Y_ts(N) = tot_ts(N, "Y")
Z_ts(N) = tot_ts(N, "Z")

"""
    X(N), Y(N), Z(N)
Total X, Y, Z operators for a system of size `N`.
"""
X(N) = Operator(X_ts(N))
Y(N) = Operator(Y_ts(N))
Z(N) = Operator(Z_ts(N))



"""
    XYZ_ts(N, deltay, deltaz)

1D XYZ model
``H= \\sum_{j=1}^N ( X_j X_{j+1} + \\Delta_Y Y_j Y_{j+1} + \\Delta_Z Z_j Z_{j+1})``
Returns a `OperatorTS1D` object.
"""
function XYZ_ts(N, deltay, deltaz)
    H = Operator(N)
    H += "X", 1, "X", 2
    H += deltay, "Y", 1, "Y", 2
    H += deltaz, "Z", 1, "Z", 2
    return OperatorTS1D(H, full=false)
end

"""
    XXZ_ts(N, delta)

1D XXZ model
``H= \\sum_{j=1}^N ( X_j X_{j+1} + Y_j Y_{j+1} + \\Delta Z_j Z_{j+1})``
Returns a `OperatorTS1D` object.
"""
XXZ_ts(N, delta) = XYZ_ts(N, 1, delta)



"""
    XY_ts(N, delta)

1D XY model
``H= \\sum_{j=1}^N ( X_j X_{j+1} +  \\Delta Y_j Y_{j+1} )``
Returns a `OperatorTS1D` object.
"""
XY_ts(N, delta) = XYZ_ts(N, delta, 0)


"""
    XX_ts(N)

1D XX model
``H= \\sum_{j=1}^N ( X_j X_{j+1} +  Y_j Y_{j+1} )``
Returns a `OperatorTS1D` object.
"""
XX_ts(N) = XYZ_ts(N, 1, 0)


"""
    XXX_ts(N)

1D XXX model
``H= \\sum_{j=1}^N ( X_j X_{j+1} + Y_j Y_{j+1} + Z_j Z_{j+1})``
Returns a `OperatorTS1D` object.
"""
XXX_ts(N) = XYZ_ts(N, 1, 1)


end
