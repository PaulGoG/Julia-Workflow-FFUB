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

function Y_A(dy, A)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.A_H):maximum(dy.A_H)
        Suma_Y = sum(dy.Y[dy.A_H .== i])
        Suma_σ = sqrt(sum(dy.σY[dy.A_H .== i].^2))
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

    f = 200/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f
    return Y
end

function Y_Z(dy, Z)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.Z_H):maximum(dy.Z_H)
        Suma_Y = sum(dy.Y[dy.Z_H .== i])
        Suma_σ = sqrt(sum(dy.σY[dy.Z_H .== i].^2))
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

    f = 200/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f    
    return Y
end

function Y_N(dy, A, Z)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.A_H):maximum(dy.A_H)
        for j in minimum(dy.Z_H[dy.A_H .== i]):maximum(dy.Z_H[dy.A_H .== i])
            Suma_Y = sum(dy.Y[dy.A_H .- dy.Z_H .== i - j])
            Suma_σ = sqrt(sum(dy.σY[dy.A_H .- dy.Z_H .== i - j].^2))
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

    f = 200/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f    
    return Y
end

function Y_TKE(dy)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.TKE):maximum(dy.TKE)
        Suma_Y = sum(dy.Y[dy.TKE .== i])
        Suma_σ = sqrt(sum(dy.σY[dy.TKE .== i].^2))
        push!(Y.x, i)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
    end

    f = 100/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f    
    return Y
end

function TKE_A(dy)
    TKE = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.A_H):maximum(dy.A_H)
        Numarator = 0
        Numitor = sum(dy.Y[dy.A_H .== i])
        Suma_σ = 0
        for j in minimum(dy.TKE):maximum(dy.TKE)
            Y_A_TKE = sum(dy.Y[(dy.A_H .== i) .& (dy.TKE .== j)])
            Numarator += j * Y_A_TKE
        end
        for j in minimum(dy.TKE):maximum(dy.TKE)
            σY_A_TKE = sqrt(sum(dy.σY[(dy.A_H .== i) .& (dy.TKE .== j)].^2))
            Suma_σ += (j*Numitor - Numarator)^2 * σY_A_TKE^2
        end
        push!(TKE.x, i)
        push!(TKE.y, Numarator/Numitor)
        push!(TKE.σ, sqrt(Suma_σ/Numitor^4))
    end  
    return TKE
end

function KE_A(tke_A, A)
    KE = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.A_H):maximum(dy.A_H)
        TKE_A = tke_A.y[tke_A.x .== i][1]
        σ_TKE_A = tke_A.σ[tke_A.x .== i][1]
        KE_H = TKE_A * (A - i)/A
        KE_L = TKE_A * i/A

        push!(KE.x, i)
        push!(KE.y, KE_H)
        push!(KE.σ, KE_H * σ_TKE_A/TKE_A)
        push!(KE.x, A - i)
        push!(KE.y, KE_L)
        push!(KE.σ, KE_L * σ_TKE_A/TKE_A)
    end
    if iseven(A)
        index = findfirst(x -> x == Int(A/2), KE.x)
        deleteat!(KE.x, index)
        deleteat!(KE.y, index)
        deleteat!(KE.σ, index)
    end
    return KE
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

function Medie_distributie_Y(distributie, index_min, index_max)
    Numarator = 0
    Numitor = 0
    for i in index_min:index_max
        Numarator = Numarator + distributie.y[i]*distributie.x[i]
        Numitor = Numitor + distributie.y[i]
    end
    return round(Numarator/Numitor, digits = 3)
end

function Medie_distributie(distributie, Y, index_min, index_max)
    Numarator = 0
    Numitor = 0
    for i in index_min:index_max
        Numarator += distributie.y[i] * Y.y[Y.x .== distributie.x[i]][1]
        Numitor += Y.y[Y.x .== distributie.x[i]][1]
    end
    return round(Numarator/Numitor, digits = 3)
end
# Aici se opreste partea de calcul a programului

# Constructia reprezentarilor grafice
function Grafic_scatter(distributie, titlu, axa_x, axa_y, scalare)
    plt = scatter(
        distributie.x, 
        distributie.y, 
        yerr = distributie.σ, 
        xlims = (minimum(distributie.x), maximum(distributie.x)),
        ylims = (minimum(distributie.y), maximum(distributie.y)*scalare),
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

y_A = Sortare_distributie(Y_A(dy, A₀))
y_Z = Sortare_distributie(Y_Z(dy, Z₀))
y_N = Sortare_distributie(Y_N(dy, A₀, Z₀))
y_TKE = Sortare_distributie(Y_TKE(dy))
tke_A = Sortare_distributie(TKE_A(dy))
ke_A = Sortare_distributie(KE_A(tke_A, A₀))

Plot_Y_A = Grafic_scatter(y_A, "Y(A)", "A", "Y %", 1.1);
Plot_Y_A = Grafic_unire_linie(y_A, Plot_Y_A);
mid_index = Int((length(y_A.x) + 1 )/2);
A_L_mediu = Medie_distributie_Y(y_A, 1, mid_index);
A_H_mediu = Medie_distributie_Y(y_A, mid_index, length(y_A.x));
Plot_Y_A = Grafic_linie_medie(Plot_Y_A, A_L_mediu);
Plot_Y_A = Grafic_linie_medie(Plot_Y_A, A_H_mediu);
Plot_Y_A = Grafic_textbox_medie(y_A.x[mid_index], maximum(y_A.y), Plot_Y_A, "A_H", A_H_mediu, 1e-3, "");
Plot_Y_A = Grafic_textbox_medie(y_A.x[mid_index], maximum(y_A.y)*0.95, Plot_Y_A, "A_L", A_L_mediu, 1e-3, "");
Grafic_afisare(Plot_Y_A, "Y(A)");

Plot_Y_Z = Grafic_scatter(y_Z, "Y(Z)", "Z", "Y %", 1.1);
Plot_Y_Z = Grafic_unire_linie(y_Z, Plot_Y_Z);
mid_index = Int((length(y_Z.x) + 1 )/2)
δₑₒ = round((sum(y_Z.y[iseven.(y_Z.x)]) - sum(y_Z.y[isodd.(y_Z.x)]))/sum(y_Z.y), digits = 3)
Plot_Y_Z = Grafic_textbox(y_Z.x[mid_index], maximum(y_Z.y), Plot_Y_Z, "\\delta_{eo}", δₑₒ, 1e-3, "");
Grafic_afisare(Plot_Y_Z, "Y(Z)");

Plot_Y_N = Grafic_scatter(y_N, "Y(N)", "N", "Y %", 1.1);
Plot_Y_N = Grafic_unire_linie(y_N, Plot_Y_N);
Grafic_afisare(Plot_Y_N, "Y(N)");

Plot_Y_TKE = Grafic_scatter(y_TKE, "Y(TKE)", "TKE [MeV]", "Y %", 1.1);
Plot_Y_TKE = Grafic_unire_linie(y_TKE, Plot_Y_TKE);
TKE_mediu = Medie_distributie_Y(y_TKE, firstindex(y_TKE.x), lastindex(y_TKE.x));
mid_index = Int((length(y_TKE.x) + 1 )/2);
Plot_Y_TKE = Grafic_linie_medie(Plot_Y_TKE, TKE_mediu);
Plot_Y_TKE = Grafic_textbox_medie(y_TKE.x[mid_index]*0.95, maximum(y_TKE.y), Plot_Y_TKE, "TKE", TKE_mediu, 1e-3, "");
Grafic_afisare(Plot_Y_TKE, "Y(TKE)");

Plot_TKE_A = Grafic_scatter(tke_A, "TKE(A)", latexstring("\$\\mathrm{A_H}\$"), "TKE [MeV]", 1.05);
TKE_A_mediu = Medie_distributie(tke_A, y_A, firstindex(tke_A.x), lastindex(tke_A.x));
Plot_TKE_A = Grafic_unire_linie(tke_A, Plot_TKE_A);
mid_index = Int((length(tke_A.x) + 1 )/2);
Grafic_textbox_medie(tke_A.x[mid_index]*1.05, maximum(tke_A.y), Plot_TKE_A, "TKE", TKE_A_mediu, 1e-3, "");
Grafic_afisare(Plot_TKE_A, "TKE(A)")

Plot_KE_A = Grafic_scatter(ke_A, "KE(A)", "A", "KE [MeV]", 1.05);
Grafic_afisare(Plot_KE_A, "KE(A)")