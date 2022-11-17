using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul pentru distributii uzuale in fisiunea U5

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

# Citire fisiere de date
df = CSV.File("Data_files/AUDI95.csv"; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σD"]) |> DataFrame
dy = CSV.File("Data_files/U5YAZTKE.csv"; delim=' ', ignorerepeated=true, header=["A_H", "Z_H", "TKE", "Y", "σY"]) |> DataFrame

struct distributie_unidym
    x
    y
    σ
end

function Y_A(dy, A, f)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.A_H):maximum(dy.A_H)
        Suma_Y = sum(dy.Y[dy.A_H .== i]) * f
        Suma_σ = sqrt(sum(dy.σY[dy.A_H .== i].^2)) * f
        push!(Y.x, i)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
        push!(Y.x, A - i)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
    end
    if iseven(A)
        index = findfirst(x -> x == Int(A/2), Y.x)
        deleteat!(Y.x, index)
        deleteat!(Y.y, index)
        deleteat!(Y.σ, index)
    end
    return Y
end

function Y_Z(dy, Z, f)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.Z_H):maximum(dy.Z_H)
        Suma_Y = sum(dy.Y[dy.Z_H .== i]) * f
        Suma_σ = sqrt(sum(dy.σY[dy.Z_H .== i].^2)) * f
        push!(Y.x, i)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
        push!(Y.x, Z - i)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
    end
    index_true = unique(i -> Y.x[i], eachindex(Y.x))
    index_delete = setdiff(eachindex(Y.x), index_true)
    deleteat!(Y.x, index_delete)
    deleteat!(Y.y, index_delete)
    deleteat!(Y.σ, index_delete)
    return Y
end

function Y_N(dy, A, Z, f)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.A_H):maximum(dy.A_H)
        for j in minimum(dy.Z_H[dy.A_H .== i]):maximum(dy.Z_H[dy.A_H .== i])
            Suma_Y = sum(dy.Y[dy.A_H .- dy.Z_H .== i - j]) * f
            Suma_σ = sqrt(sum(dy.σY[dy.A_H .- dy.Z_H .== i - j].^2)) * f
            push!(Y.x, i - j)
            push!(Y.y, Suma_Y)
            push!(Y.σ, Suma_σ)
            push!(Y.x, A - Z - i + j)
            push!(Y.y, Suma_Y)
            push!(Y.σ, Suma_σ)
        end
    end
    index_true = unique(i -> Y.x[i], eachindex(Y.x))
    index_delete = setdiff(eachindex(Y.x), index_true)
    deleteat!(Y.x, index_delete)
    deleteat!(Y.y, index_delete)
    deleteat!(Y.σ, index_delete)
    return Y
end
# Aici se opreste partea de calcul a programului

# Constructia reprezentarilor grafice
function Grafic_scatter(distributie, titlu)
    plt = scatter(
        distributie.x, 
        distributie.y, 
        yerr = distributie.σ, 
        xlims = (minimum(distributie.x), maximum(distributie.x)),
        ylims = (minimum(distributie.y), maximum(distributie.y)*1.1),
        xlabel = "X axis", 
        ylabel = "Y axis", 
        framestyle = :box,
        legend = :false,
        title = "$titlu",
        minorgrid = :true,
        size = (900, 900)
    )
    return plt
end

function Grafic_adauga_linie(distributie, plt, titlu, val_medie, sigma_medie)
    x = sort(distributie.x)
    y = [distributie.y[distributie.x .== i][1] for i in minimum(x):maximum(x)]
    σ = [distributie.σ[distributie.x .== i][1] for i in minimum(x):maximum(x)]
    plot!(plt, 
        x, 
        y,
        ribbon = σ
    )
    plt = annotate!(maximum(x)*0.925, maximum(y)*0.99, latexstring("\$<\\mathrm{Q}> =\$ $val_medie \$\\pm\$ $sigma_medie MeV"))
    plt = vline!(x[y .== val_medie], ls = :dashdot, label = "")
    display(plt)
    savefig(plt, "Grafice/$(titlu)_T2.png")
end

# Apelarea functiilor definite pentru executia programului
A = 236
Z = 92
f = 100/sum(dy.Y)

y_A = Y_A(dy, A, f)
y_Z = Y_Z(dy, Z, f)
y_N = Y_N(dy, A, Z, f)

Grafic_adauga_linie(y_A, Grafic_scatter(y_A, "Y(A)"), "Y(A)", 0, 0)
Grafic_adauga_linie(y_Z, Grafic_scatter(y_Z, "Y(Z)"), "Y(Z)", 0, 0)
Grafic_adauga_linie(y_N, Grafic_scatter(y_N, "Y(N)"), "Y(N)", 0, 0)

#Efect Even-Odd
(sum(y_Z.y[iseven.(y_Z.x)]) - sum(y_Z.y[isodd.(y_Z.x)]))/sum(y_Z.y)