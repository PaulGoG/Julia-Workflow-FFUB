using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul pentru distributii uzuale in fisiune

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

struct distributie_bidym
    x_1
    x_2
    y
    σ
end

function Y_A(dy, A)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for A_H in minimum(dy.A_H):maximum(dy.A_H)
        # Y(A) = Σ_(Z, TKE) Y(A, Z, TKE)
        # σY(A) = sqrt[Σ_(Z, TKE) σY(A, Z, TKE)^2]
        Suma_Y = sum(dy.Y[dy.A_H .== A_H])
        Suma_σ = sqrt(sum(dy.σY[dy.A_H .== A_H].^2))
        push!(Y.x, A_H)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
        if A - A_H != A_H
            push!(Y.x, A - A_H)
            push!(Y.y, Suma_Y)
            push!(Y.σ, Suma_σ)
        end
    end
    # Normarea distributiei
    f = 200/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f
    return Y
end

function Y_Z(dy, Z)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for A_H in minimum(dy.Z_H):maximum(dy.Z_H)
        # Y(Z) = Σ_(A, TKE) Y(A, Z, TKE)
        # σY(Z) = sqrt[Σ_(A, TKE) σY(A, Z, TKE)^2]
        Suma_Y = sum(dy.Y[dy.Z_H .== A_H])
        Suma_σ = sqrt(sum(dy.σY[dy.Z_H .== A_H].^2))
        push!(Y.x, A_H)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
        push!(Y.x, Z - A_H)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
    end
    # Stergerea datelor duplicate
    index_true = unique(i -> Y.x[i], eachindex(Y.x))
    index_delete = setdiff(eachindex(Y.x), index_true)
    deleteat!(Y.x, index_delete)
    deleteat!(Y.y, index_delete)
    deleteat!(Y.σ, index_delete)
    # Normarea distributiei
    f = 200/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f    
    return Y
end

function Y_N(dy, A, Z)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for A_H in minimum(dy.A_H):maximum(dy.A_H)
        for Z_H in minimum(dy.Z_H[dy.A_H .== A_H]):maximum(dy.Z_H[dy.A_H .== A_H])
            Suma_Y = sum(dy.Y[dy.A_H .- dy.Z_H .== A_H - Z_H])
            Suma_σ = sqrt(sum(dy.σY[dy.A_H .- dy.Z_H .== A_H - Z_H].^2))
            push!(Y.x, A_H - Z_H)
            push!(Y.y, Suma_Y)
            push!(Y.σ, Suma_σ)
            push!(Y.x, A - Z - A_H + Z_H)
            push!(Y.y, Suma_Y)
            push!(Y.σ, Suma_σ)
        end
    end
    # Stergerea datelor duplicate
    index_true = unique(i -> Y.x[i], eachindex(Y.x))
    index_delete = setdiff(eachindex(Y.x), index_true)
    deleteat!(Y.x, index_delete)
    deleteat!(Y.y, index_delete)
    deleteat!(Y.σ, index_delete)
    # Normarea distributiei
    f = 200/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f    
    return Y
end

function Y_TKE(dy)
    Y = distributie_unidym(Int[],Float64[], Float64[])
    for TKE in minimum(dy.TKE):maximum(dy.TKE)
        # Y(TKE) = Σ_(A, Z) Y(A, Z, TKE)
        # σY(TKE) = sqrt[Σ_(A, Z) σY(A, Z, TKE)^2]
        Suma_Y = sum(dy.Y[dy.TKE .== TKE])
        Suma_σ = sqrt(sum(dy.σY[dy.TKE .== TKE].^2))
        push!(Y.x, TKE)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
    end
    # Normarea distributiei
    f = 100/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f    
    return Y
end

function TKE_A(dy)
    tke = distributie_unidym(Int[],Float64[], Float64[])
    for A_H in minimum(dy.A_H):maximum(dy.A_H)
        Numarator = 0
        Numitor = sum(dy.Y[dy.A_H .== A_H])
        Suma_σ² = 0
        # TKE(A) = Σ_(TKE) TKE * Y(A, TKE)/Σ_(TKE) Y(A, TKE)
        # σTKE(A) = [1/Σ_(TKE) Y(A, TKE)] * sqrt[Σ_(TKE) (TKE - TKE(A)) * σY(A, TKE)^2]
        for TKE in minimum(dy.TKE[dy.A_H .== A_H]):maximum(dy.TKE[dy.A_H .== A_H])
            Y_A_TKE = sum(dy.Y[(dy.A_H .== A_H) .& (dy.TKE .== TKE)])
            Numarator += TKE * Y_A_TKE
        end
        tke_A = Numarator/Numitor
        for TKE in minimum(dy.TKE[dy.A_H .== A_H]):maximum(dy.TKE[dy.A_H .== A_H])
            σY_A_TKE = sqrt(sum(dy.σY[(dy.A_H .== A_H) .& (dy.TKE .== TKE)].^2))
            Suma_σ² += (TKE - tke_A)^2 * σY_A_TKE^2
        end
        push!(tke.x, A_H)
        push!(tke.y, tke_A)
        push!(tke.σ, sqrt(Suma_σ²)/Numitor)
    end  
    return tke
end

function KE_A(tke_A, A)
    KE = distributie_unidym(Int[],Float64[], Float64[])
    for A_H in minimum(dy.A_H):maximum(dy.A_H)
        TKE_A = tke_A.y[tke_A.x .== A_H][1]
        σTKE_A = tke_A.σ[tke_A.x .== A_H][1]
        # KEₗₕ(A) = Aₕₗ/A₀ * TKE(Aₕ)
        # σKEₗₕ(A) = Aₕₗ/A₀ * σTKE(Aₕ)
        KE_H = TKE_A * (A - A_H)/A
        KE_L = TKE_A * A_H/A
        push!(KE.x, A_H)
        push!(KE.y, KE_H)
        push!(KE.σ, KE_H * σTKE_A/TKE_A)
        if A - A_H != A_H
            push!(KE.x, A - A_H)
            push!(KE.y, KE_L)
            push!(KE.σ, KE_L * σTKE_A/TKE_A)
        end
    end
    return KE
end

function Q_A_Z(A, Z, df, limInfA_H, limSupA_H)
    Q = distributie_bidym(Int[], Int[], Float64[], Float64[])
    D = df.D[(df.A .== A) .& (df.Z .== Z)][1]
    σ_D = df.σD[(df.A .== A) .& (df.Z .== Z)][1]
    for A_H in limInfA_H:limSupA_H
        for Z_H in minimum(df.Z[df.A .== A_H]):maximum(df.Z[df.A .== A_H])
            if isassigned(df.D[(df.A .== A_H) .& (df.Z .== Z_H)], 1) && isassigned(df.D[(df.A .== A - A_H) .& (df.Z .== Z - Z_H)], 1)
                D_H = df.D[(df.A .== A_H) .& (df.Z .== Z_H)][1]
                σ_D_H = df.σD[(df.A .== A_H) .& (df.Z .== Z_H)][1]
                D_L = df.D[(df.A .== A - A_H) .& (df.Z .== Z - Z_H)][1]
                σ_D_L = df.σD[(df.A .== A - A_H) .& (df.Z .== Z - Z_H)][1]
                q = D - (D_H + D_L)
                # Energiile sunt salvate în MeV
                push!(Q.y, q *1e-3)
                push!(Q.σ, sqrt(σ_D^2 + σ_D_H^2 + σ_D_L^2) *1e-3)
                push!(Q.x_1, A_H)
                push!(Q.x_2, Z_H)
            end
        end
    end 
    return Q
end

function Q_A(q_A_Z, y_Z)
    Q = distributie_unidym(Int[], Float64[], Float64[])
    for A_H in minimum(q_A_Z.x_1):maximum(q_A_Z.x_1)
        Numarator = 0
        Numitor = 0
        Suma_σ² = 0
        # Q(A) = Σ_(Z) Q(A, Z) * Y(Z)/Σ_(Z) Y(Z)
        if isassigned(q_A_Z.x_2[q_A_Z.x_1 .== A_H], 1)
            for Z_H = minimum(q_A_Z.x_2[q_A_Z.x_1 .== A_H]):maximum(q_A_Z.x_2[q_A_Z.x_1 .== A_H])
                if isassigned(y_Z.y[y_Z.x .== Z_H], 1)
                    Numarator += q_A_Z.y[(q_A_Z.x_1 .== A_H) .& (q_A_Z.x_2 .== Z_H)][1] * y_Z.y[y_Z.x .== Z_H][1]
                    Numitor += y_Z.y[y_Z.x .== Z_H][1]
                    Suma_σ² += y_Z.y[y_Z.x .== Z_H][1]^2 * q_A_Z.σ[(q_A_Z.x_1 .== A_H) .& (q_A_Z.x_2 .== Z_H)][1]^2
                end
            end
            if Numitor != 0 # Elementele (A, Z) pentru care a existat și Q(A, Z) și Y(Z)
                q_A = Numarator/Numitor
                for Z_H = minimum(q_A_Z.x_2[q_A_Z.x_1 .== A_H]):maximum(q_A_Z.x_2[q_A_Z.x_1 .== A_H])
                    if isassigned(y_Z.y[y_Z.x .== Z_H], 1)
                        Suma_σ² += y_Z.σ[y_Z.x .== Z_H][1]^2 * (q_A_Z.y[(q_A_Z.x_1 .== A_H) .& (q_A_Z.x_2 .== Z_H)][1] - q_A)^2
                    end
                end
                push!(Q.y, q_A)
                push!(Q.σ, sqrt(Suma_σ²)/Numitor)
                push!(Q.x, A_H)
            end
        end
    end
    return Q
end

# Calculul energiei de separare a particulei (A_part, Z_part) din nucleul (A, Z)
function Energie_separare(A_part, Z_part, A, Z, df)
    # Verificarea existentei radionuclizilor folositi in libraria de date
    if isassigned(df.D[(df.A .== A_part) .& (df.Z .== Z_part)], 1) && isassigned(df.D[(df.A .== A) .& (df.Z .== Z)], 1) && isassigned(df.D[(df.A .== A - A_part) .& (df.Z .== Z - Z_part)], 1)
        D_part = df.D[(df.A .== A_part) .& (df.Z .== Z_part)][1]
        σ_part = df.σD[(df.A .== A_part) .& (df.Z .== Z_part)][1]
    
        D = df.D[(df.A .== A - A_part) .& (df.Z .== Z - Z_part)][1]
        σᴰ = df.σD[(df.A .== A - A_part) .& (df.Z .== Z - Z_part)][1]    
    
        S = (D + D_part - df.D[(df.A .== A) .& (df.Z .== Z)][1])*1e-3
        σˢ = sqrt(σ_part^2 + σᴰ^2 + (df.σD[(df.A .== A) .& (df.Z .== Z)][1])^2)*1e-3

        return [S, σˢ]
    else 
        return [NaN, NaN]
    end 
end

function TXE_A(q_A, tke_A, df, A, Z)
    TXE = distributie_unidym(Int[],Float64[], Float64[])
    Sₙ = Energie_separare(1, 0, A, Z, df)
    for A_H in minimum(q_A.x):maximum(q_A.x)
        if isassigned(tke_A.y[tke_A.x .== A_H], 1)
            # TXE(A) = Q(A) - TKE(A) + Sₙ
            push!(TXE.x, A_H)
            push!(TXE.y, q_A.y[q_A.x .== A_H][1] + Sₙ[1] - tke_A.y[tke_A.x .== A_H][1])
            push!(TXE.σ, sqrt(q_A.σ[q_A.x .== A_H][1]^2 + Sₙ[2]^2 + tke_A.σ[tke_A.x .== A_H][1]^2))
        end
    end  
    return TXE
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
    # Valoare medie & incertitudine pentru distributiile de yield
    Numarator = 0
    Numitor = 0
    Suma_σ² = 0
    for i in index_min:index_max
        Numarator += distributie.y[i]*distributie.x[i]
        Numitor += distributie.y[i]
    end
    Media_distributiei = Numarator/Numitor
    for i in index_min:index_max
        Suma_σ² += (distributie.x[i] - Media_distributiei)^2 * distributie.σ[i]^2
    end    
    return [round(Media_distributiei, digits = 3), round(sqrt(Suma_σ²)/Numitor, digits = 5)]
end

function Medie_distributie(distributie, Y, index_min, index_max)
    # Valoare medie & incertitudine pentru distributii mediate folosind o distributie de yield
    Numarator = 0
    Numitor = 0
    Suma_σ² = 0
    for i in index_min:index_max
        if isassigned(Y.y[Y.x .== distributie.x[i]], 1)
            Numarator += distributie.y[i] * Y.y[Y.x .== distributie.x[i]][1]
            Numitor += Y.y[Y.x .== distributie.x[i]][1]
        end
    end
    if Numitor != 0 # Elementele pentru care a existat și Xᵢ și Yᵢ la același indice
        Media_distributiei = Numarator/Numitor
        for i in index_min:index_max
            if isassigned(Y.y[Y.x .== distributie.x[i]], 1)
                Suma_σ² += Y.y[Y.x .== distributie.x[i]][1]^2 * distributie.σ[i]^2
                Suma_σ² += (distributie.x[i] - Media_distributiei)^2 * Y.σ[Y.x .== distributie.x[i]][1]^2
            end
        end    
        return [round(Media_distributiei, digits = 3), round(sqrt(Suma_σ²)/Numitor, digits = 5)]    
    else 
        return [NaN, NaN]
    end
end
# Aici se opreste partea de calcul a programului

# Constructia reprezentarilor grafice
function Grafic_scatter(distributie, titlu, axa_x, axa_y, scalare_inf, scalare_sup)
    plt = scatter(
        distributie.x, 
        distributie.y, 
        yerr = distributie.σ, 
        xlims = (minimum(distributie.x), maximum(distributie.x)),
        ylims = (minimum(distributie.y)*scalare_inf, maximum(distributie.y)*scalare_sup),
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
    fillalpha = .3,
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
function Grafic_linie_medie_vertical(plt, distributie_med)
    vline!(plt, 
    [distributie_med], 
    ls = :dashdot, 
    label = ""
    )
    return plt
end
function Grafic_linie_medie_orizontal(plt, distributie_med)
    hline!(plt, 
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
limInfA_H = minimum(dy.A_H)
limSupA_H = maximum(dy.A_H)

y_A = Sortare_distributie(Y_A(dy, A₀))
y_Z = Sortare_distributie(Y_Z(dy, Z₀))
y_N = Sortare_distributie(Y_N(dy, A₀, Z₀))
y_TKE = Sortare_distributie(Y_TKE(dy))
tke_A = Sortare_distributie(TKE_A(dy))
ke_A = Sortare_distributie(KE_A(tke_A, A₀))
q_A = Sortare_distributie(Q_A(Q_A_Z(A₀, Z₀, df, limInfA_H, limSupA_H), y_Z))
txe_A = Sortare_distributie(TXE_A(q_A, tke_A, df, A₀, Z₀))

Plot_Y_A = Grafic_scatter(y_A, "Y(A)", "A", "Y %", 1, 1.1);
Plot_Y_A = Grafic_unire_linie(y_A, Plot_Y_A);
mid_index = Int((length(y_A.x) + 1 )/2);
A_L_mediu = Medie_distributie_Y(y_A, 1, mid_index);
A_H_mediu = Medie_distributie_Y(y_A, mid_index, length(y_A.x));
Plot_Y_A = Grafic_linie_medie_vertical(Plot_Y_A, A_L_mediu[1]);
Plot_Y_A = Grafic_linie_medie_vertical(Plot_Y_A, A_H_mediu[1]);
Plot_Y_A = Grafic_textbox_medie(y_A.x[mid_index], maximum(y_A.y), Plot_Y_A, "A_H", A_H_mediu[1], A_H_mediu[2], "");
Plot_Y_A = Grafic_textbox_medie(y_A.x[mid_index], maximum(y_A.y)*0.95, Plot_Y_A, "A_L", A_L_mediu[1], A_L_mediu[2], "");
Grafic_afisare(Plot_Y_A, "Y(A)");

Plot_Y_Z = Grafic_scatter(y_Z, "Y(Z)", "Z", "Y %", 1, 1.1);
Plot_Y_Z = Grafic_unire_linie(y_Z, Plot_Y_Z);
mid_index = Int((length(y_Z.x) + 1 )/2)
δₑₒ = (sum(y_Z.y[iseven.(y_Z.x)]) - sum(y_Z.y[isodd.(y_Z.x)]))/sum(y_Z.y)
σ_δₑₒ = round((1/sum(y_Z.y)) * sqrt((1 + δₑₒ)^2 * sum(y_Z.σ .^2) + 2*δₑₒ*(sum(y_Z.σ[isodd.(y_Z.x)] .^2) - sum(y_Z.σ[iseven.(y_Z.x)].^2))), digits = 7)
δₑₒ = round(δₑₒ, digits = 5)
Plot_Y_Z = Grafic_textbox(y_Z.x[mid_index], maximum(y_Z.y), Plot_Y_Z, "\\delta_{eo}", δₑₒ, σ_δₑₒ, "");
Grafic_afisare(Plot_Y_Z, "Y(Z)");

Plot_Y_N = Grafic_scatter(y_N, "Y(N)", "N", "Y %", 1, 1.1);
Plot_Y_N = Grafic_unire_linie(y_N, Plot_Y_N);
Grafic_afisare(Plot_Y_N, "Y(N)");

Plot_Y_TKE = Grafic_scatter(y_TKE, "Y(TKE)", "TKE [MeV]", "Y %", 1, 1.1);
Plot_Y_TKE = Grafic_unire_linie(y_TKE, Plot_Y_TKE);
TKE_mediu = Medie_distributie_Y(y_TKE, firstindex(y_TKE.x), lastindex(y_TKE.x));
mid_index = Int((length(y_TKE.x) + 1 )/2);
Plot_Y_TKE = Grafic_linie_medie_vertical(Plot_Y_TKE, TKE_mediu[1]);
Plot_Y_TKE = Grafic_textbox_medie(y_TKE.x[mid_index]*0.94, maximum(y_TKE.y), Plot_Y_TKE, "TKE", TKE_mediu[1], TKE_mediu[2], "MeV");
Grafic_afisare(Plot_Y_TKE, "Y(TKE)");

Plot_TKE_A = Grafic_scatter(tke_A, "TKE(A)", latexstring("\$\\mathrm{A_H}\$"), "TKE [MeV]", 1, 1.05);
TKE_A_mediu = Medie_distributie(tke_A, y_A, firstindex(tke_A.x), lastindex(tke_A.x));
Plot_TKE_A = Grafic_unire_linie(tke_A, Plot_TKE_A);
mid_index = Int((length(tke_A.x) + 1 )/2);
Plot_TKE_A = Grafic_textbox_medie(tke_A.x[mid_index]*1.05, maximum(tke_A.y), Plot_TKE_A, "TKE", TKE_A_mediu[1], TKE_A_mediu[2], "MeV");
Grafic_afisare(Plot_TKE_A, "TKE(A)")

Plot_KE_A = Grafic_scatter(ke_A, "KE(A)", "A", "KE [MeV]", 1, 1.05);
Plot_KE_A = Grafic_unire_linie(ke_A, Plot_KE_A);
Grafic_afisare(Plot_KE_A, "KE(A)")

Plot_Q_A = Grafic_scatter(q_A, "Q(A)", latexstring("\$\\mathrm{A_H}\$"), "Q [MeV]", 0.98, 1.02);
Plot_Q_A = Grafic_unire_linie(q_A, Plot_Q_A);
Q_A_Mediu = Medie_distributie(q_A, y_A, firstindex(q_A.x), lastindex(q_A.x));
mid_index = Int((length(q_A.x) + 1)/2);
Plot_Q_A = Grafic_textbox_medie(q_A.x[mid_index], maximum(q_A.y), Plot_Q_A, "Q", Q_A_Mediu[1], Q_A_Mediu[2], "MeV");
Plot_Q_A = Grafic_linie_medie_orizontal(Plot_Q_A, Q_A_Mediu[1]);
Grafic_afisare(Plot_Q_A, "Q(A)")

Plot_TXE_A = Grafic_scatter(txe_A, "TXE(A)", latexstring("\$\\mathrm{A_H}\$"), "TXE [MeV]", 0.9, 1.1);
Plot_TXE_A = Grafic_unire_linie(txe_A, Plot_TXE_A);
mid_index = Int((length(txe_A.x) + 1)/2);
TXE_A_Mediu = Medie_distributie(txe_A, y_A, firstindex(txe_A.x), lastindex(txe_A.x));
Plot_TXE_A = Grafic_textbox_medie(txe_A.x[mid_index], maximum(txe_A.y), Plot_TXE_A, "TXE", TXE_A_Mediu[1], TXE_A_Mediu[2], "MeV");
Plot_TXE_A = Grafic_linie_medie_orizontal(Plot_TXE_A, TXE_A_Mediu[1]);
Grafic_afisare(Plot_TXE_A, "TXE(A)")