module ModelsOpSum

using ITensors, ITensorMPS


function ising(N, hx, hz)
    H = OpSum()
    for j in 1:N
        H += 4 , "Sz", j, "Sz", mod1(j + 1, N)
        H += 2*hx, "Sx", j
        H += 2*hz, "Sz", j
    end
    return H
end


end
