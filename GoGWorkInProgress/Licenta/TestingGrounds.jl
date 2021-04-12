using Plots
using CSV
using DataFrames

include("Constante.jl")
include("RawDataRead.jl")
include("TestDependency.jl")
#include("Helpers.jl")

S = 0
for k in 1:16
    F_k = first(Freq.F_k[(Freq.Zona_k .== k)])
    F_ki = Freq.F_ki[(Freq.Zona_k .== k)]
    global S = S + F_k * sum(F_ki)
end

Suprafata = "Agricol"
H_1(Suprafata)
