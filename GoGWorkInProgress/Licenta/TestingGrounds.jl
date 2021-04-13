using Plots
using CSV
using DataFrames

include("Constante.jl")
include("RawDataRead.jl")
include("TestingDependency.jl")

#=
Plotare Cladiri
x=y = collect(-1000.0:1:1000.0)
unghi = collect(0:0.01:2*pi)
Cercx = 1000 *cos.(unghi)
Cercy = 1000 *sin.(unghi)
scatter(Cladiri.x, Cladiri.y,Cladiri.z)
=#
