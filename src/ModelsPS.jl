
module ModelsPS

using PauliStrings


function ising_ts(N, hx, hz)
    H = Operator(N)
    H += "Z", 1, "Z", 2
    H += hx, "X", 1
    H += hz, "Z", 1
    return OperatorTS1D(H, full=false)
end

ising(N, hx, hz) = Operator(ising_ts(N, hx, hz))

end
