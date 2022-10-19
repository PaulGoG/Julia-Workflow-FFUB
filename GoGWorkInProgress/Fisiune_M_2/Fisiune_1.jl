using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul al energiei eliberate la fisiune 

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

struct Q
    Z
    A
    Q
end

struct Q_mediat
    A
    Q
end

function Q_A_Z(librarie, A, Z, limInfA_H, limSupA_H)
    q = Q(Int[], Int[], Float64[])

    # Citim fisierul de tip CSV
    # Z|A|Simbol|D(KeV)|σᴰ(KeV) -> forma tabelului
    df = CSV.File(librarie; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) |> DataFrame
    D = df.D[(df.A .== A) .& (df.Z .== Z)][1]

    # Parcurgem radionuclizii din baza de date relevanti pentru fragmentari
    for i in limInfA_H:limSupA_H
        A_H = i
        Z_UCD = Z*A_H/A
        Z_p = floor(Z_UCD - 0.5)

        for j in (Z_p - 1):(Z_p + 1)
            Z_H = j
            if isassigned(df.D[(df.A .== A_H) .& (df.Z .== Z_H)], 1) && isassigned(df.D[(df.A .== A - A_H) .& (df.Z .== Z - Z_H)], 1)
                D_H = df.D[(df.A .== A_H) .& (df.Z .== Z_H)][1]
                D_L = df.D[(df.A .== A - A_H) .& (df.Z .== Z - Z_H)][1]
                push!(q.Q, (D - (D_H + D_L))*1e-3)
                push!(q.Z, Z_H)
                push!(q.A, A_H)
            end
        end
    end 
    return q
end

function p_A_Z(Z, Z_p)
    return 1/(sqrt(2*pi) * 0.6) * exp(-(Z - Z_p)^2 /(2*0.6^2))
end
function Q_A(Q, A, Z, limInfA_H, limSupA_H)
    Q_med = Q_mediat(Int[], Float64[])

    for i in limInfA_H:limSupA_H
        A_H = i
        Mediere_sus = 0
        Mediere_jos = 0
        Z_UCD = Z*A_H/A
        Z_p = Z_UCD - 0.5
        for j = 1:length(Q.Z[Q.A .== A_H])
            Mediere_sus = Mediere_sus + Q.Q[Q.A .== A_H][j] * p_A_Z(Q.Z[Q.A .== A_H][j], Z_p)
            Mediere_jos = Mediere_jos + p_A_Z(Q.Z[Q.A .== A_H][j], Z_p)
        end
        push!(Q_med.Q, Mediere_sus/Mediere_jos)
        push!(Q_med.A, A_H)
    end
    return Q_med
end
# Aici se opreste partea de calcul a programului

# Constructia reprezentarilor grafice
# Scatter simplu
function Grafic_simplu(Q, librarie)
    plt = scatter(
        Q.A, 
        Q.Q,  
        marker = :circle,
        markersize = 5, 
        markerstrokewidth = 0,
        xlabel = L"\mathrm{A_H}", 
        ylabel = latexstring("\$Q\$  [MeV]"), 
        framestyle = :box,
        legend = :false,
        title = "Energia eliberată la fisiune, $(librarie[begin:end-4])",
        minorgrid = :true,
        mc = :red
    )
    display(plt)
    #savefig(plt, "Grafice\\Q_A_$(librarie[begin:end-4]).png")
end

# Apelarea functiilor definite pentru executia programului
audi95 = "AUDI95.csv"
audi21 = "AUDI2021.csv"
A = 236
Z = 92
limInfA_H = 118
limSupA_H = 160

Q_95_A_Z = Q_A_Z(audi95, A, Z, limInfA_H, limSupA_H)
Grafic_simplu(Q_95_A_Z, audi95)
Q_95_A = Q_A(Q_95_A_Z, A, Z, limInfA_H, limSupA_H)
Grafic_simplu(Q_95_A, audi95)