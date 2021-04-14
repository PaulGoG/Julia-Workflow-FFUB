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
plot(x,x-> H_final(x, Surface))
#hbtrans = x -> filter!(y->y!==nothing, [Δh_b_tranzitie(x[i]) for i in 1:length(x)])
#xprim = collect(0:1:length(hbtrans(x))-1)
#plot(xprim, hbtrans(x) , legend = :outertopright);
#plot(x,x-> Δh_m_tranzitie(x));
#hline!(x,x-> Δh_b_final()[1]);
#hline!(x,x-> Δh_b_final()[2]);
#hline!(x,x-> Δh_b_final()[3])
#hline!(x,x-> Δh_m_final()[1]);
#hline!(x,x-> Δh_m_final()[2]);
#hline!(x,x-> Δh_m_final()[3])