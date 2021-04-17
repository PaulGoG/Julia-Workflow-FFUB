#=
Simulam dispersia atmosferica a gazului de tritiu care iese la
cosul unui reactor nuclear de tip CANDU
!Important! -> x>0 !!!
=#

using Plots; plotlyjs()
using Trapz
using DataFrames
using CSV

t_R = 60*60*365*17 # Timpul de emisie in secunde
t_zile = t_R/86400
t_spalare = t_R * (15/365) # Consideram ca avem 15 zile ploioase intr-un an

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
dim_transversal = 10000
dim_vertical = 0
step = 1000

x = collect(-(dim_transversal/2):step:(dim_transversal/2))
y = collect(-(dim_transversal/2):step:(dim_transversal/2))
z = collect(0:1:dim_vertical)

if t_R <= 3600
    χ_Q = χ_Q_Durata_Scurta(x, y, z, Pasquill, Suprafata, Tip_Suprafata, t_R)[:,:,1]
    χ = χ_Durata_Scurta(χ_Q, x, y, z, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)[:,:,1]
    ω = ω_Scurt(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_R, t_spalare)
elseif t_R <= 86400
    χ_Q = χ_Q_Durata_Prelungita(x, y, Pasquill, Suprafata, Tip_Suprafata)
    χ = χ_Durata_Prelungita(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
    ω = ω_Prelungit(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_R, t_spalare)
else
    χ_Q = χ_Q_Durata_Lunga(x, y, Suprafata, Tip_Suprafata)
    χ = χ_Durata_Lunga(χ_Q, x, y, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
    ω = ω_Lung(χ_Q, x, y, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
end
