
using ITensors, ITensorMPS
using QuantumDynamics
using LinearAlgebra
import PyPlot as plt

N = 6
dt = 0.1
tmax = 20.0
hx = 0.5
hz = 0.0

# time evolution with dense matrices
H = Matrix(ModelsPS.ising(N, hx, hz))
O = Matrix(ModelsPS.Z(N))
traceO0(Ot) = tr(Matrix{Complex}(Ot * O))
times, res = heisenberg(H, O, [traceO0], dt, tmax)
res = real(res[1, :])
plt.plot(times, res, label="dense")



# time evolution with TEBD
sites = siteinds("S=1/2", N)
H = ModelsTEBD.ising(N, hx, hz, sites)
O = MPO(ModelsOpSum.Z(N), sites)
traceO0(Ot) = tr(apply(Ot, swapprime(dag(O), 0 => 1)))
times, res = heisenberg(H, O, [traceO0], dt, tmax)
res = real(res[1, :])
plt.plot(times, res, label="TEBD")
plt.legend()
plt.show()
