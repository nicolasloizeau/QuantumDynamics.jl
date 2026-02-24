
using LinearAlgebra
using ProgressBars
using SparseArrays



function shrodinger(H::AbstractMatrix, ψ::AbstractVector, observables, dt::Real, tmax::Real)
    U = exp(-im * H * dt)
    times = 0:dt:tmax
    res = zeros(Complex{Float64}, length(times), length(observables))
    for (i, t) in ProgressBar(enumerate(times))
        res[i, :] = [obs(ψ) for obs in observables]
        ψ = U * ψ
    end
    return times, transpose(res)
end


function rk4(H::AbstractMatrix, ψ::AbstractVector, dt::Real)
    k1 = -im * H * ψ
    k2 = -im * H * (ψ + dt/2 * k1)
    k3 = -im * H * (ψ + dt/2 * k2)
    k4 = -im * H * (ψ + dt * k3)
    return dt/6 * (k1 + 2*k2 + 2*k3 + k4)
end

function shrodinger_nh(H::AbstractMatrix, ψ::AbstractVector, observables, dt::Real, tmax::Real)
    times = 0:dt:tmax
    res = zeros(Complex{Float64}, length(times), length(observables))
    for (i, t) in ProgressBar(enumerate(times))
        res[i, :] = [obs(ψ) for obs in observables]
        ψ += rk4(H, ψ, dt)
    end
    return times, transpose(res)
end
