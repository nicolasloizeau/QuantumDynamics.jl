

using LinearAlgebra
using ProgressBars
using SparseArrays
using Expokit


"""
    heisenberg(H::AbstractMatrix, O::Matrix, observables, dt::Real, tmax::Real)

Heisenberg evolution with dense matrices.
Evolve the operator `O` in the Heisenberg picture under the Hamiltonian `H` for a time `tmax` with time step `dt`.
`observables` is a list of functions, which are applied to the operator `O` at each time step.
"""
function heisenberg(H::Matrix, O::Matrix, observables, dt::Real, tmax::Real)
    U = exp(-im * H * dt)
    times = 0:dt:tmax
    res = zeros(Complex{Float64}, length(times), length(observables))
    for (i, t) in ProgressBar(enumerate(times))
        res[i, :] = [obs(O) for obs in observables]
        O = U * O * U'
    end
    return times, transpose(res)
end


function heisenberg(H::Matrix, O::Matrix, t::Real)
    U = exp(-im * H * t)
    return U * O * U'
end



function rk4_vec(L::AbstractMatrix, Ovec::AbstractVector, dt::Real)
    k1 = L * Ovec
    k2 = L * (Ovec + dt/2 * k1)
    k3 = L * (Ovec + dt/2 * k2)
    k4 = L * (Ovec + dt * k3)
    return Ovec + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
end

function heisenberg(H::SparseMatrixCSC, O::AbstractMatrix, observables, dt::Real, tmax::Real)
    Ovec = Vector{ComplexF64}(vec(O))
    n = norm(Ovec)
    one = sparse(I, size(H, 1), size(H, 1))
    L = -im * (kron(one, H) - kron(H', one))
    times = 0:dt:tmax
    res = zeros(Complex{Float64}, length(times), length(observables))
    for (i, t) in ProgressBar(enumerate(times))
        res[i, :] = [obs(reshape(Ovec, size(O))) for obs in observables]
        # Ovec = rk4_vec(L, Ovec, dt)
        Ovec = expmv!( Ovec, dt, L, Ovec; m=5)

        # Ovec *= n / norm(Ovec)
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
    for (i, t) in ProgressBar(enumerate(times))
        res[i, :] = [obs(O) for obs in observables]
        O = apply(Hgates, O; cutoff=cutoff, apply_dag=true, maxdim=bd)
        O = truncate(O; maxdim=bd, cutoff=cutoff)
    end
    return times, transpose(res)
end

function heisenberg(H::AbstractOperator, O::AbstractOperator, observables, dt::Real, tmax::Real; M=2^16, noise=0)
    times = 0:dt:tmax
    res = zeros(Complex{Float64}, length(times), length(observables))
    for (i, t) in ProgressBar(enumerate(times))
        res[i, :] = [obs(O) for obs in observables]
        O = rk4(H, O, dt; heisenberg=true, M=M) # perform one step of rk4, keep only M strings
        O = add_noise(O, noise * dt) # add depolarizing noise
        O = trim(O, M) # keep the M strings with the largest weight
    end
    return times, transpose(res)
end


function heisenberg(H::AbstractOperator, O::AbstractOperator, dt::Real, tmax::Real; M=2^16, noise=0)
    times = 0:dt:tmax
    for (i, t) in ProgressBar(enumerate(times))
        O = rk4(H, O, dt; heisenberg=true, M=M) # perform one step of rk4, keep only M strings
        O = add_noise(O, noise * dt) # add depolarizing noise
        O = trim(O, M) # keep the M strings with the largest weight
    end
    return O
end
