#=
Simulam dispersia atmosferica a gazului de tritiu care iese la
cosul unui reactor nuclear de tip CANDU
!Important! -> x>0 !!!
=#

using Plots; plotlyjs()
using Trapz
using DataFrames
using CSV

const t_R = 86400*365 # Timpul de emisie in secunde
const C = 2 # Intervine la calculul Σ-urilor, trebuie sa maximizeze dilutia

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
dim_transversal = 100000
dim_vertical = 0
step = 1000

x = collect(-(dim_transversal/2):step:(dim_transversal/2))
y = collect(-(dim_transversal/2):step:(dim_transversal/2))
z = collect(0:1:dim_vertical)

if t_R <= 3600
    χ_Q = χ_Q_Durata_Scurta(x, y, z, Pasquill, Suprafata, Tip_Suprafata)
    χ = χ_Durata_Scurta(χ_Q, x, y, z, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
    ω = ω_Scurt(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
elseif t_R <= 86400
    χ_Q = χ_Q_Durata_Prelungita(x, y, Pasquill, Suprafata, Tip_Suprafata)
    χ = χ_Durata_Prelungita(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
    ω = ω_Prelungit(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
else
    χ_Q = χ_Q_χ_Durata_Lunga(x, y, Suprafata, Tip_Suprafata)
    χ = χ_Durata_Lunga(χ_Q, x, y, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
    ω = ω_Lung(χ_Q, x, y, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
end

Reprezentari(x, y, χ_Q)
Reprezentari(x, y, χ)
Reprezentari(x, y, ω)
Reprezentari(x, y, ω .* Resuspensie())