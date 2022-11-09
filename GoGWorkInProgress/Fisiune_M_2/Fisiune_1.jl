using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul al energiei eliberate la fisiune folosind 3 Z/A in jurul Z_UCD

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

struct distributie_bidym
    x_1
    x_2
    y
    σ
end

# Functie de calcul pentru Q(A, Z)
function Q_A_Z(librarie, A, Z, limInfA_H, limSupA_H)
    Q = distributie_bidym(Int[], Int[], Float64[], Float64[])

    # Citim fisierul de tip CSV
    # Z|A|Simbol|D(KeV)|σᴰ(KeV) -> forma tabelului
    df = CSV.File(librarie; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) |> DataFrame
    D = df.D[(df.A .== A) .& (df.Z .== Z)][1]
    σ_D = df.σ[(df.A .== A) .& (df.Z .== Z)][1]

    for i in limInfA_H:limSupA_H
        A_H = i
        Z_UCD = Z*A_H/A
        Z_p = floor(Z_UCD - 0.5)
        for j in (Z_p - 1):(Z_p + 1) # Fragmentarile 3 Z/A
            Z_H = j
            if isassigned(df.D[(df.A .== A_H) .& (df.Z .== Z_H)], 1) && isassigned(df.D[(df.A .== A - A_H) .& (df.Z .== Z - Z_H)], 1)
                D_H = df.D[(df.A .== A_H) .& (df.Z .== Z_H)][1]
                σ_D_H = df.σ[(df.A .== A_H) .& (df.Z .== Z_H)][1]
                D_L = df.D[(df.A .== A - A_H) .& (df.Z .== Z - Z_H)][1]
                σ_D_L = df.σ[(df.A .== A - A_H) .& (df.Z .== Z - Z_H)][1]
                push!(Q.y, (D - (D_H + D_L)) *1e-3)
                push!(Q.σ, sqrt(σ_D^2 + σ_D_H^2 + σ_D_L^2) *1e-3)
                push!(Q.x_1, A_H)
                push!(Q.x_2, Z_H)
            end
        end
    end 
    return Q
end

# Distributia izobara de sarcina cu rms(A) ≈ 0.6 & ΔZₚ ≈ 0.5
function p_A_Z(Z, Z_p)
    return 1/(sqrt(2*π) * 0.6) * exp(-(Z - Z_p)^2 /(2*0.6^2))
end

function Q_A(Q, A, Z, limInfA_H, limSupA_H)
    Q_med = distributie_bidym(Int[], Int[], Float64[], Float64[])

    for i in limInfA_H:limSupA_H
        A_H = i
        Numarator = 0
        Numitor = 0
        Sigma_temp² = 0
        Z_UCD = Z*A_H/A
        Z_p = Z_UCD - 0.5
        for j = 1:length(Q.x_2[Q.x_1 .== A_H])
            Numarator = Numarator + Q.y[Q.x_1 .== A_H][j] * p_A_Z(Q.x_2[Q.x_1 .== A_H][j], Z_p)
            Numitor = Numitor + p_A_Z(Q.x_2[Q.x_1 .== A_H][j], Z_p)
            Sigma_temp² = Sigma_temp² + (p_A_Z(Q.x_2[Q.x_1 .== A_H][j], Z_p) * Q.σ[Q.x_1 .== A_H][j])^2
        end
        push!(Q_med.y, Numarator/Numitor)
        push!(Q_med.σ, sqrt(Sigma_temp²)/Numitor)
        push!(Q_med.x_1, A_H)
    end
    return Q_med
end
# Aici se opreste partea de calcul a programului

# Constructia reprezentarilor grafice
function Grafic_scatter(Q, eticheta)
    plt = scatter(
        Q.x_1, 
        Q.y,  
        yerr = Q.σ, 
        markersize = 4, 
        markerstrokewidth = 1,
        xlims = (minimum(Q.x_1) - 1, maximum(Q.x_1) + 1),
        xlabel = L"\mathrm{A_H}", 
        ylabel = latexstring("\$\\mathrm{Q}\$  [MeV]"), 
        framestyle = :box,
        label = "$eticheta",
        title = "Energia eliberată la fisiune",
        minorgrid = :true,
        msc = :red,
        size = (800, 600)
    )
    return plt
    #savefig(plt, "Grafice\\NumeGrafic.png")
end
function Grafic_suprapunere(Q, plt, eticheta)
    scatter!(plt, 
    Q.x_1, 
    Q.y,
    yerr = Q.σ,
    label = "$eticheta",
    markersize = 4, 
    markerstrokewidth = 1,
    msc = :black
)
display(plt)
end

# Apelarea functiilor definite pentru executia programului
audi95 = "AUDI95.csv"
A = 236
Z = 92
limInfA_H = 118
limSupA_H = 160

QAZ = Q_A_Z(audi95, A, Z, limInfA_H, limSupA_H)
QA = Q_A(QAZ, A, Z, limInfA_H, limSupA_H)
Q_mediu = sum(QA.y)/length(QA.y)
σ_Q_mediu = 1/length(QA.y) * sqrt(sum(QA.σ .^2))
Grafic_suprapunere(QA, Grafic_scatter(QAZ, "Q(A, Z)"), "Q(A)")