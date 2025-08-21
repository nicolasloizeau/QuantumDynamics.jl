

using LinearAlgebra
using ProgressBars



"""
    heisenberg(H::Matrix, O::Matrix, observables, dt::Real, tmax::Real)

Heisenberg evolution with dense matrices.
Evolve the operator `O` in the Heisenberg picture under the Hamiltonian `H` for a time `tmax` with time step `dt`.
`observables` is a list of functions, which are applied to the operator `O` at each time step.
"""
function heisenberg(H::Matrix, O::Matrix, observables, dt::Real, tmax::Real)
    U = exp(-im * H * dt)
    times = 0:dt:tmax
    res = zeros(Complex{Float64}, length(times), length(observables))
    for (i,t) in ProgressBar(enumerate(times))
        res[i,:] = [obs(O) for obs in observables]
        O = U * O * U'
    end
    return times, transpose(res)
end

function togates(terms, step)
    gates = [exp(-im * step / 2 * term) for term in terms]
    append!(gates, reverse(gates))
end

"""
    heisenberg(H::Vector{ITensor}, O::MPO, observables, dt::Real, tmax::Real; bd=100, cutoff=1E-8)

Heisenberg evolution with MPO (TEBD).
Evolve the operator `O` in the Heisenberg picture under the Hamiltonian `H` for a time `tmax` with time step `dt`.
`observables` is a list of functions, which are applied to the operator `O` at each time step.
`H` is a list of ITensors representing the Hamiltonian terms.
"""
function heisenberg(H::Vector{ITensor}, O::MPO, observables, dt::Real, tmax::Real; bd=100, cutoff=1E-8)
    Hgates = togates(H, dt)
    times = 0:dt:tmax
    res = zeros(Complex{Float64}, length(times), length(observables))
    for (i,t) in ProgressBar(enumerate(times))
        res[i,:] = [obs(O) for obs in observables]
        O = apply(Hgates, O; cutoff=cutoff, apply_dag=true, maxdim=bd)
        O = truncate(O; maxdim=bd, cutoff=cutoff)
    end
    return times, transpose(res)
end
