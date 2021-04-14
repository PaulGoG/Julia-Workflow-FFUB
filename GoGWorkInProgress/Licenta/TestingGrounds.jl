using Plots; gr()
using CSV
using DataFrames

include("Constante.jl")
include("CitireDate.jl")
include("TestingDependency.jl")

#=
Plotare Cladiri
x=y = collect(-1000.0:1:1000.0)
unghi = collect(0:0.01:2*pi)
Cercx = 1000 *cos.(unghi)
Cercy = 1000 *sin.(unghi)
scatter(Cladiri.x, Cladiri.y,Cladiri.z)
=#

x = collect(0:1:1000.0)
Surface = unique(T_4.Tip_Suprafata)[rand((1,2,3))] # For testing purposes
Pasquill = "D"
plot(x,x-> H_final(x, Pasquill, Surface))

