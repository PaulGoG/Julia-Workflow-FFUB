using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul pentru energia eliberata la fisiune folosind 3 Z/A fragmente in jurul Zₚ(A)

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

struct distributie_bidym
    x_1
    x_2
    y
    σ
end
#####
# Functie de calcul pentru Q(A,Z)
function Q_A_Z(librarie, A, Z, limInfA_H, limSupA_H)
    Q = distributie_bidym(Int[], Int[], Float64[], Float64[])
    # Citim fisierul de tip CSV
    # Z|A|Simbol|D(KeV)|σᴰ(KeV) -> forma tabelului
    df = CSV.read(librarie, DataFrame; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"])
    D = df.D[(df.A .== A) .& (df.Z .== Z)][1]
    σ_D = df.σ[(df.A .== A) .& (df.Z .== Z)][1]
    for A_H in limInfA_H:limSupA_H
        Z_UCD = Z*A_H/A
        Z_p = Z_UCD - 0.5
        Z_mid = round(Z_p)
        for Z_H in (Z_mid - 1):(Z_mid + 1) # Fragmentarile 3 Z/A
            # Verificam ca exista radionuclizii in libraria de date
            A_L = A - A_H
            Z_L = Z - Z_H
            if isassigned(df.D[(df.A .== A_H) .& (df.Z .== Z_H)], 1) && isassigned(df.D[(df.A .== A_L) .& (df.Z .== Z_L)], 1)
                D_H = df.D[(df.A .== A_H) .& (df.Z .== Z_H)][1]
                σ_D_H = df.σ[(df.A .== A_H) .& (df.Z .== Z_H)][1]
                D_L = df.D[(df.A .== A_L) .& (df.Z .== Z_L)][1]
                σ_D_L = df.σ[(df.A .== A_L) .& (df.Z .== Z_L)][1]
                # Valorile sunt calculate in MeVi
                push!(Q.y, (D - (D_H + D_L)) *1e-3)
                push!(Q.σ, sqrt(σ_D^2 + σ_D_H^2 + σ_D_L^2) *1e-3)
                push!(Q.x_1, A_H)
                push!(Q.x_2, Z_H)
            end
        end
    end 
    return Q
end
# Distributia izobara de sarcina cu  σ = rms(A) ≈ 0.6 & ΔZₚ ≈ 0.5 (Gaussiana)
function p_A_Z(Z, Z_p)
    return 1/(sqrt(2*π) * 0.6) * exp(-(Z - Z_p)^2 /(2* 0.6^2))
end
# Functia de mediere a Q(A,Z) pe distributia p(A,Z) considerand 3 Z/A fragmente
function Q_A(q_A_Z, A, Z)
    q_A = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in minimum(q_A_Z.x_1):maximum(q_A_Z.x_1)
        Numarator = 0
        Numitor = 0
        Sigma_temp² = 0
        Z_UCD = Z*A_H/A
        Z_p = Z_UCD - 0.5
        # Q(A) = Σ_Z Q(A, Z)*p(A, Z)/Σ_Z p(A, Z)
        # σQ(A) = sqrt[Σ_Z σQ(A, Z)^2 *p(A, Z)^2]/Σ_Z p(A, Z)
        for i = 1:length(q_A_Z.x_2[q_A_Z.x_1 .== A_H])
            P_A_Z = p_A_Z(q_A_Z.x_2[q_A_Z.x_1 .== A_H][i], Z_p)
            Numarator += q_A_Z.y[q_A_Z.x_1 .== A_H][i] * P_A_Z
            Numitor += P_A_Z
            Sigma_temp² += (P_A_Z * q_A_Z.σ[q_A_Z.x_1 .== A_H][i])^2
        end
        if Numitor != 0
            push!(q_A.y, Numarator/Numitor)
            push!(q_A.σ, sqrt(Sigma_temp²)/Numitor)
            push!(q_A.x_1, A_H)   
        end
    end
    return q_A
end
# Aici se opreste partea de calcul a programului
#####
# Constructia reprezentarilor grafice
function Grafic_scatter(Q, eticheta, culoare)
    plt = scatter(
        Q.x_1, 
        Q.y,  
        yerr = Q.σ, 
        xlabel = L"\mathrm{A_H}", 
        ylabel = "Q [MeV]", 
        xlims = (minimum(Q.x_1)-0.5, maximum(Q.x_1)+0.5),
        framestyle = :box,
        label = eticheta,
        title = latexstring("Energia eliberată la fisiune calculată folosind 3 Z/A în jurul \$\\mathrm{Z_{p}(A)}\$ fragmente"),
        minorgrid = :true,
        size = (1000, 950),
        dpi = 600,
        color = culoare,
        markersize = 6
    )
    return plt
end
function Grafic_scatter(Q, plt, eticheta, culoare)
    scatter!(plt, 
        Q.x_1, 
        Q.y,  
        yerr = Q.σ,
        label = eticheta,
        color = culoare,
        markersize = 6
    )
    return plt
end
function Grafic_unire_linie(Q, plt, culoare)
    plot!(plt, 
    Q.x_1, 
    Q.y,
    ribbon = Q.σ,
    #fillalpha = .3,
    label = "",
    color = culoare
    )
    return plt
end
function Grafic_afisare(plt)
    display(plt)
    savefig(plt, "Grafice/T1_Q_value.png")
end
#####
# Apelarea functiilor definite pentru executia programului
librarie = "Data_files/Defecte_masa/AUDI2021.csv"
A₀ = 236;
Z₀ = 92;
limInfA_H = 118;
limSupA_H = 160;

q_A_Z = Q_A_Z(librarie, A₀, Z₀, limInfA_H, limSupA_H);
q_A = Q_A(q_A_Z, A₀, Z₀);

Plot = Grafic_scatter(q_A_Z, "Q(A,Z)", :navy);
Plot = Grafic_scatter(q_A, Plot, "Q(A)", :red);
Plot = Grafic_unire_linie(q_A, Plot, :red);
Grafic_afisare(Plot);