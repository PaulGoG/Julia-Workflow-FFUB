#=
Simulam dispersia atmosferica a gazului de tritiu care iese la
cosul unui reactor nuclear de tip CANDU
!Important! -> x>0 !!!
=#

using Plots; plotlyjs()
using Trapz
using DataFrames
using CSV

const t_R = 3000 # Timpul de emisie in secunde
const C = 2 # Intervine la calculul Σ-urilor, normal trebuie sa maximizeze Dilutia

include("Constante.jl")
include("CitireDate.jl")
include("Helpers.jl")
include("Calcul_dilutie.jl")
include("Vectorize.jl")
include("ReprezentariGrafice.jl")

Pasquill = "D"
Suprafata = "Padure_Oras"
Tip_Suprafata = "Padure_Urban"
Tip_Aversa = "Ploaie"
Debit = 1.0
dim_transversal = 1000
dim_vertical = 0

x = collect(-(dim_transversal/2):1:(dim_transversal/2))
y = collect(-(dim_transversal/2):1:(dim_transversal/2))
z = collect(0:1:dim_vertical)

Aux = Durata_Scurta(x, y, z, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
χ_Q = Aux[1][:,:,1]
χ = Aux[2][:,:,1]
ω = Aux[3]
K = Aux[4]