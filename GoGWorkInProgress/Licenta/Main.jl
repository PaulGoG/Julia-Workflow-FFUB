#=
Simulam dispersia atmosferica a gazului de tritiu care iese la
cosul unui reactor nuclear de tip CANDU
!Important! -> x>0 !!!
=#

using Plots; plotlyjs()
using Trapz
using DataFrames
using CSV

t_R = 8640*365/2 # Timpul de emisie in secunde
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
Debit = 3.0
dim_transversal = 10000
dim_vertical = 0
step = 100

x = collect(-(dim_transversal/2):step:(dim_transversal/2))
y = collect(-(dim_transversal/2):step:(dim_transversal/2))
z = collect(0:1:dim_vertical)

if t_R <= 3600
    χ_Q = χ_Q_Durata_Scurta(x, y, z, Pasquill, Suprafata, Tip_Suprafata, t_R)[:,:,1]
    χ = χ_Durata_Scurta(χ_Q, x, y, z, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)[:,:,1]
    ω = ω_Scurt(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_R, t_spalare)

    Reprezinta_Suprafata(x, y, χ_Q, x[Int(dim_transversal/(step*2) + 1)], "χ/Q (s/m^3)", "χ_Q")
    Reprezinta_Gradient(x, y, χ_Q, x[Int(dim_transversal/(step*2) + 1)], "χ/Q (s/m^3)", "χ_Q")
    Reprezinta_Contur(x, y, χ_Q, x[Int(dim_transversal/(step*2) + 1)], "χ/Q (s/m^3)", "χ_Q")

    Reprezinta_Suprafata(x, y, χ, x[Int(dim_transversal/(step*2) + 1)], "χ (Bq s/m^3)", "χ")
    Reprezinta_Gradient(x, y, χ, x[Int(dim_transversal/(step*2) + 1)], "χ (Bq s/m^3)", "χ")
    Reprezinta_Contur(x, y, χ, x[Int(dim_transversal/(step*2) + 1)], "χ (Bq s/m^3)", "χ")

    Reprezinta_Suprafata(x, y, ω, x[Int(dim_transversal/(step*2) + 1)], "ω (Bq/m^2)", "ω")
    Reprezinta_Gradient(x, y, ω, x[Int(dim_transversal/(step*2) + 1)], "ω (Bq/m^2)", "ω")
    Reprezinta_Contur(x, y, ω, x[Int(dim_transversal/(step*2) + 1)], "ω (Bq/m^2)", "ω")

    K = Coeficient_Resuspensie(t_zile)
    Resuspensie = ω .* K
    Reprezinta_Suprafata(x, y, Resuspensie, x[Int(dim_transversal/(step*2) + 1)], "Resuspensie (Bq/m^3)", "Resuspensie")
    Reprezinta_Gradient(x, y, Resuspensie, x[Int(dim_transversal/(step*2) + 1)], "Resuspensie (Bq/m^3)", "Resuspensie")
    Reprezinta_Contur(x, y, Resuspensie, x[Int(dim_transversal/(step*2) + 1)], "Resuspensie (Bq/m^3)", "Resuspensie")
elseif t_R <= 86400
    χ_Q = χ_Q_Durata_Prelungita(x, y, Pasquill, Suprafata, Tip_Suprafata)
    χ = χ_Durata_Prelungita(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
    ω = ω_Prelungit(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_R, t_spalare)

    Reprezinta_Suprafata(x, y, χ_Q, x[Int(dim_transversal/(step*2) + 1)], "χ/Q (s/m^3)", "χ_Q")
    Reprezinta_Gradient(x, y, χ_Q, x[Int(dim_transversal/(step*2) + 1)], "χ/Q (s/m^3)", "χ_Q")
    Reprezinta_Contur(x, y, χ_Q, x[Int(dim_transversal/(step*2) + 1)], "χ/Q (s/m^3)", "χ_Q")

    Reprezinta_Suprafata(x, y, χ, x[Int(dim_transversal/(step*2) + 1)], "χ (Bq s/m^3)", "χ")
    Reprezinta_Gradient(x, y, χ, x[Int(dim_transversal/(step*2) + 1)], "χ (Bq s/m^3)", "χ")
    Reprezinta_Contur(x, y, χ, x[Int(dim_transversal/(step*2) + 1)], "χ (Bq s/m^3)", "χ")

    Reprezinta_Suprafata(x, y, ω, x[Int(dim_transversal/(step*2) + 1)], "ω (Bq/m^2)", "ω")
    Reprezinta_Gradient(x, y, ω, x[Int(dim_transversal/(step*2) + 1)], "ω (Bq/m^2)", "ω")
    Reprezinta_Contur(x, y, ω, x[Int(dim_transversal/(step*2) + 1)], "ω (Bq/m^2)", "ω")

    K = Coeficient_Resuspensie(t_zile)
    Resuspensie = ω .* K
    Reprezinta_Suprafata(x, y, Resuspensie, x[Int(dim_transversal/(step*2) + 1)], "Resuspensie (Bq/m^3)", "Resuspensie")
    Reprezinta_Gradient(x, y, Resuspensie, x[Int(dim_transversal/(step*2) + 1)], "Resuspensie (Bq/m^3)", "Resuspensie")
    Reprezinta_Contur(x, y, Resuspensie, x[Int(dim_transversal/(step*2) + 1)], "Resuspensie (Bq/m^3)", "Resuspensie")
else
    χ_Q = χ_Q_Durata_Lunga(x, y, Suprafata, Tip_Suprafata)
    χ = χ_Durata_Lunga(χ_Q, x, y, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
    ω = ω_Lung(χ_Q, x, y, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)

    Reprezinta_Suprafata(x, y, χ_Q, x[1], "χ/Q (s/m^3)", "χ_Q")
    Reprezinta_Gradient(x, y, χ_Q, x[1], "χ/Q (s/m^3)", "χ_Q")
    Reprezinta_Contur(x, y, χ_Q, x[1], "χ/Q (s/m^3)", "χ_Q")

    Reprezinta_Suprafata(x, y, χ, x[1], "χ (Bq s/m^3)", "χ")
    Reprezinta_Gradient(x, y, χ, x[1], "χ (Bq s/m^3)", "χ")
    Reprezinta_Contur(x, y, χ, x[1], "χ (Bq s/m^3)", "χ")

    Reprezinta_Suprafata(x, y, ω, x[1], "ω (Bq/m^2)", "ω")
    Reprezinta_Gradient(x, y, ω, x[1], "ω (Bq/m^2)", "ω")
    Reprezinta_Contur(x, y, ω, x[1], "ω (Bq/m^2)", "ω")

    K = Coeficient_Resuspensie(t_zile)
    Resuspensie = ω .* K
    Reprezinta_Suprafata(x, y, Resuspensie, x[1], "Resuspensie (Bq/m^3)", "Resuspensie")
    Reprezinta_Gradient(x, y, Resuspensie, x[1], "Resuspensie (Bq/m^3)", "Resuspensie")
    Reprezinta_Contur(x, y, Resuspensie, x[1], "Resuspensie (Bq/m^3)", "Resuspensie")
end
