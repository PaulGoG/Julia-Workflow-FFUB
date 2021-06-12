#=
Fisierul principal al programului;
Kernelul trebuie oprit la fiecare modificare a timpului 
de emisie t_R deoarece prin intermediul lui se defineste 
activitatea totala a poluantului emis Q_0 in Constante.jl
=#

using Plots; plotlyjs()
using Trapz
using DataFrames
using CSV

t_R = 3600*24*365*5 # Timpul de emisie in secunde
t_zile = t_R/86400
# In medie avem 450.6 ore de precipitatii pe an
t_spalare = t_R * (450.6/(365*24))

include("Constante.jl")
include("CitireDate.jl")
include("Helpers.jl")
include("Calcul_dilutie.jl")
include("Vectorize.jl")
include("ReprezentariGrafice.jl")

#= 
Pentru emisiile cu t < 24 ore folosim clasa de stabilitate atmosferica 
neutra deoarece are frecventa de aparitie de aproape 50%
=#
Pasquill = "D"; 

# Suprafetele depind de specificul zonei
Suprafata = "Agricol"; # Tabel_4
Tip_Suprafata = "Rural"; # Tabel_2

# Aversa dominanta este ploaia (> 70% din timp)
Tip_Aversa = "Ploaie"
Debit = 0.5

# Dimensiunile zonei analizate in m
dim_transversal = 60000
step = 150 # Pasul de evaluare
dim_vertical = 0

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