using Plots;
using CSV
using DataFrames

include("Constante.jl")
include("CitireDate.jl")
include("Helpers.jl")

#=
Plotare Cladiri
x=y = collect(-1000.0:1:1000.0)
unghi = collect(0:0.01:2*pi)
Cercx = 1000 *cos.(unghi)
Cercy = 1000 *sin.(unghi)
scatter(Cladiri.x, Cladiri.y,Cladiri.z)
=#

#=
x = collect(0:1:1000.0)
Surface = unique(T_4.Tip_Suprafata)[rand((1,2,3))] # For testing purposes
Pasquill = "D"
plot(x,x-> H_final(x, Pasquill, Surface))
=#

#=
p = dilutie_instantanee(100000000, 15, 0, Pasquill, Suprafata, Tip_Suprafata)
C = 0.5
C = 1
C = 2
x = collect(1:1:100)
sy = [σ_y(x[i], Pasquill) for i in 1:length(x)]
sz = [σ_z(x[i], Pasquill, Tip_Suprafata) for i in 1:length(x)]
Sz = [Σ_z(x[i], Pasquill, Suprafata, Tip_Suprafata) for i in 1:length(x)]
Sy = [Σ_y(x[i], Pasquill, Suprafata) for i in 1:length(x)]
plot(x,sy);
plot!(x,sz);
plot!(x,Sz);
plot!(x,Sy)
=#