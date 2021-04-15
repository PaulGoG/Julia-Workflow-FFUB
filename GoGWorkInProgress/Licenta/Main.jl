#=
Simulam dispersia atmosferica a gazului de tritiu care iese la
cosul unui reactor nuclear de tip CANDU
!Important! -> x>0 !!!
=#

using Plots; plotlyjs()
using Trapz
using DataFrames
using CSV

t_R = 86400E2 # Timpul de emisie in secunde

include("Constante.jl")
include("CitireDate.jl")
include("Helpers.jl")
include("Calcul_dilutie.jl")
include("ReprezentariGrafice.jl")


Pasquill = "C"
Suprafata = "Padure_Oras"
Tip_Suprafata = "Padure_Urban"
dim_transversal = 1000
dim_vertical = 0

x = collect(1:1:dim_transversal+1)
y = collect(-(dim_transversal/2):1:(dim_transversal/2))
z = collect(0:1:dim_vertical)


χ_Q = zeros(length(x), length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            χ_Q[i,j] = dilutie_lunga_durata(x[i],y[j], Suprafata, Tip_Suprafata)
        end 
    end
heatmap(transpose(χ_Q), xaxis = "x", yaxis = "y")
contourf(transpose(χ_Q), xaxis = "x", yaxis = "y")
surface(transpose(χ_Q), xaxis = "x", yaxis = "y")