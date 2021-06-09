#=
Program de test unde verific daca algoritmul
Stie sa determine singur sectorul dintr-un arc de cerc unde
Se afla un anumit punct (x,y)
=#

using Plots; gr()
using Trapz
using DataFrames
using CSV

include("Constante.jl")
include("Helpers.jl")

θ = collect(0:0.1:2*π)
x = 100*cos.(θ)
y = 100*sin.(θ)
xprim = collect(-100:1:100)
plot(x,y, legend = false);
hline!([0]);
vline!([0]);
for l in 1:(n/2)
    if l%4 != 0
    yprim = [tan(l*θ_L)*xprim[o] for o in 1:length(xprim)]
    yprim = yprim[((yprim[:].^2 + xprim[:].^2) .< 10000)]
    xplot = collect(-length(yprim)/2:1:length(yprim)/2 -1)
    plot!(xplot,yprim)
    end
end
a = x[rand(1:1:length(x))]
b = y[rand(1:1:length(y))]
scatter!([a],[b])
Apartenenta_Sector_Cerc(a,b)