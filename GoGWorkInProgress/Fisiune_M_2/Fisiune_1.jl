using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul al energiei eliberate la fisiune folosind 3 Z/A in jurul Z_UCD

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#De schimbat numele in definitia structului la ceva mai general a.i. Sa putem folosi Q la obiectul propriu-zis
#Suprapunerea ploturilor Q(A,Z) & Q(A) cu scris pe grafice
#Propagarea erorilor

struct distributie_bidym
    x_1
    x_2
    y
    σ
end


function Q_A_Z(librarie, A, Z, limInfA_H, limSupA_H)
    Q = distributie_bidym(Int[], Int[], Float64[], Float64[])

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
function Q_A(q, A, Z, limInfA_H, limSupA_H)
    q_med = Q(Int[], Int[], Float64[], Float64[])

    for i in limInfA_H:limSupA_H
        A_H = i
        Mediere_sus = 0
        Mediere_jos = 0
        Z_UCD = Z*A_H/A
        Z_p = Z_UCD - 0.5
        for j = 1:length(q.Z[q.A .== A_H])
            Mediere_sus = Mediere_sus + q.Q[q.A .== A_H][j] * p_A_Z(q.Z[q.A .== A_H][j], Z_p)
            Mediere_jos = Mediere_jos + p_A_Z(q.Z[q.A .== A_H][j], Z_p)
        end
        push!(q_med.Q, Mediere_sus/Mediere_jos)
        push!(q_med.A, A_H)
    end
    return q_med
end
# Aici se opreste partea de calcul a programului

# Constructia reprezentarilor grafice
function Grafic_scatter(distributie)
    plt = scatter(
        Q.A, 
        Q.Q,  
        marker = :circle,
        markersize = 5, 
        markerstrokewidth = 0,
        xlabel = L"\mathrm{A_H}", 
        ylabel = latexstring("\$\\mathrm{Q}\$  [MeV]"), 
        framestyle = :box,
        legend = :false,
        title = "Energia eliberată la fisiune",
        minorgrid = :true,
        mc = :red
    )
    display(plt)
    #savefig(plt, "Grafice\\Q_A_$(librarie[begin:end-4]).png")
end

# Apelarea functiilor definite pentru executia programului
audi95 = "AUDI95.csv"
A = 236
Z = 92
limInfA_H = 118
limSupA_H = 160

Q_95_A_Z = Q_A_Z(audi95, A, Z, limInfA_H, limSupA_H)
Grafic_simplu(Q_95_A_Z)
Q_95_A = Q_A(Q_95_A_Z, A, Z, limInfA_H, limSupA_H)
Grafic_simplu(Q_95_A)