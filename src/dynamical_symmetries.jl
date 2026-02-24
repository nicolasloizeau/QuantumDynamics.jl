

function buildL(b; open=false)
    diag_values = zeros(length(b))
    L = Tridiagonal(b[1:end-1], diag_values, -b[1:end-1])
    if open
        L[end, end] = -2 * b[end]
        L[end, end-1] += b[end]
    end
    return L
end

function op_from_phi(phi, On)
    O = 0
    for n in 0:min(length(phi), length(On))-1
        O += (1im)^n * On[n+1] * phi[n+1]
    end
    return O
end

function phi_from_op(O, On)
    return [1 / (1im)^n * ps.trace_product(O, On[n+1]) for n in 0:length(On)-1]
end


function get_left_index(λr, λl; atol=1e-12, rtol=1e-8)
    left_index = similar(collect(1:length(λr)))
    for i in 1:length(λr)
        # try to find exact match (within atol/rtol), otherwise pick the nearest
        j = findfirst(k -> isapprox(λl[k], conj(λr[i]); atol=atol, rtol=rtol), 1:length(λl))
        if j === nothing
            _, j = findmin(abs.(λl .- conj(λr[i])))
        end
        left_index[i] = j
    end
    return left_index
end


"""
    dynamical_symmetries(H::AbstractOperator, seed::AbstractOperator; steps_lanczos=20, nterms_lanczos=2^20)

Compute the dynamical symmetries of the Hamiltonian `H` using the Lanczos algorithm
"""
function dynamical_symmetries(H::AbstractOperator, seed::AbstractOperator; steps_lanczos=20, nterms_lanczos=2^20)
    (bn, On) = lanczos(H, seed, steps_lanczos, nterms_lanczos; returnOn=true)
    L = buildL(bn; open=true)
    omegasr, phis = eigen(L)
    Aright = [op_from_phi(phis[:, i], On) for i in 1:length(omegasr)]
    omegasl, phis = eigen(L')
    left_index = get_left_index(omegasr, omegasl)
    Aleft = [op_from_phi(phis[:, left_index[i]], On) for i in 1:length(omegasr)]
    return Aright, Aleft, omegasr
end
