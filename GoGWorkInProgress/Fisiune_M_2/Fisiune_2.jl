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

function Y_TKE(dy, f)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.TKE):maximum(dy.TKE)
        Suma_Y = sum(dy.Y[dy.TKE .== i]) * f
        Suma_σ = sqrt(sum(dy.σY[dy.TKE .== i].^2)) * f
        push!(Y.x, i)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
    end
    return Y
end

function Sortare_distributie(distributie)
    x = sort(distributie.x)
    y = [distributie.y[distributie.x .== i][1] for i in minimum(x):maximum(x)]
    σ = [distributie.σ[distributie.x .== i][1] for i in minimum(x):maximum(x)]
    for i in eachindex(x)
        distributie.x[i] = x[i]
        distributie.y[i] = y[i]
        distributie.σ[i] = σ[i]
    end
    return distributie
end

function Medie_distributie(distributie, index_min, index_max)
    Numarator = 0
    Numitor = 0
    for i in index_min:index_max
        Numarator = Numarator + distributie.y[i]*distributie.x[i]
        Numitor = Numitor + distributie.y[i]
    end
    return round(Numarator/Numitor, digits = 3)
end
# Aici se opreste partea de calcul a programului

# Constructia reprezentarilor grafice
function Grafic_scatter(distributie, titlu, axa_x, axa_y)
    plt = scatter(
        distributie.x, 
        distributie.y, 
        yerr = distributie.σ, 
        xlims = (minimum(distributie.x), maximum(distributie.x)),
        ylims = (minimum(distributie.y), maximum(distributie.y)*1.1),
        xlabel = "$axa_x", 
        ylabel = "$axa_y", 
        framestyle = :box,
        legend = :false,
        title = "$titlu",
        minorgrid = :true,
        size = (900, 900)
    )
    return plt
end
function Grafic_unire_linie(distributie, plt)
    plot!(plt, 
    distributie.x, 
    distributie.y,
    ribbon = distributie.σ,
    label = ""
    )
    return plt
end
function Grafic_textbox_medie(x, y, plt, distributie_nume, distributie_med, distributie_med_sigma, unitate_masura)
    annotate!(plt, 
    x, 
    y, 
    latexstring("\$<\\mathrm{$distributie_nume}> =\$ $distributie_med \$\\pm\$ $distributie_med_sigma $unitate_masura")
    )
    return plt
end
function Grafic_textbox(x, y, plt, distributie_nume, distributie_val, distributie_val_sigma, unitate_masura)
    annotate!(plt, 
    x, 
    y, 
    latexstring("\$\\mathrm{$distributie_nume} =\$ $distributie_val \$\\pm\$ $distributie_val_sigma $unitate_masura")
    )
    return plt
end
function Grafic_linie_medie(plt, distributie_med)
    vline!(plt, 
    [distributie_med], 
    ls = :dashdot, 
    label = ""
    )
    return plt
end
function Grafic_afisare(plt, titlu)
    display(plt)
    #savefig(plt, "Grafice/$(titlu)_T2.png")
end

# Apelarea functiilor definite pentru executia programului
A₀ = 236
Z₀ = 92
f = 100/sum(dy.Y)

y_A = Sortare_distributie(Y_A(dy, A₀, f))
y_Z = Sortare_distributie(Y_Z(dy, Z₀, f))
y_N = Sortare_distributie(Y_N(dy, A₀, Z₀, f))
y_TKE = Sortare_distributie(Y_TKE(dy, f))

Plot_Y_A = Grafic_scatter(y_A, "Y(A)", "A", "Y %");
Plot_Y_A = Grafic_unire_linie(y_A, Plot_Y_A);
mid_index = Int((length(y_A.x) + 1 )/2);
A_L_mediu = Medie_distributie(y_A, 1, mid_index);
A_H_mediu = Medie_distributie(y_A, mid_index, length(y_A.x));
Plot_Y_A = Grafic_linie_medie(Plot_Y_A, A_L_mediu);
Plot_Y_A = Grafic_linie_medie(Plot_Y_A, A_H_mediu);
Plot_Y_A = Grafic_textbox_medie(y_A.x[mid_index], maximum(y_A.y), Plot_Y_A, "A_H", A_H_mediu, 1e-3, "");
Plot_Y_A = Grafic_textbox_medie(y_A.x[mid_index], maximum(y_A.y)*0.95, Plot_Y_A, "A_L", A_L_mediu, 1e-3, "");
Grafic_afisare(Plot_Y_A, "Y(A)");

Plot_Y_Z = Grafic_scatter(y_Z, "Y(Z)", "Z", "Y %");
Plot_Y_Z = Grafic_unire_linie(y_Z, Plot_Y_Z);
mid_index = Int((length(y_Z.x) + 1 )/2)
δₑₒ = round((sum(y_Z.y[iseven.(y_Z.x)]) - sum(y_Z.y[isodd.(y_Z.x)]))/sum(y_Z.y), digits = 3)
Plot_Y_Z = Grafic_textbox(y_Z.x[mid_index], maximum(y_Z.y), Plot_Y_Z, "\\delta_{eo}", δₑₒ, 1e-3, "");
Grafic_afisare(Plot_Y_Z, "Y(Z)");

Plot_Y_N = Grafic_scatter(y_N, "Y(N)", "N", "Y %");
Plot_Y_N = Grafic_unire_linie(y_N, Plot_Y_N);
Grafic_afisare(Plot_Y_N, "Y(N)");

Plot_Y_TKE = Grafic_scatter(y_TKE, "Y(TKE)", "TKE", "Y %");
Plot_Y_TKE = Grafic_unire_linie(y_TKE, Plot_Y_TKE);
Grafic_afisare(Plot_Y_TKE, "Y(TKE)");