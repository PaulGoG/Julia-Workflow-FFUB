using Plots
using CSV
using DataFrames
using LaTeXStrings
using QuadGK
using Trapz

# Cod de calcul pentru model de emisie cu tratare globala, 1 fragmentare cea mai probabila (Los Alamos)
# σ_c constant & P(T) parametrizat cu un triunghi dreptunghic

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

# Citire fisiere de date
df = CSV.read("Data_files/Defecte_masa/AUDI2021.csv", DataFrame; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σD"]);
dy = CSV.read("Data_files/Yield/U5YAZTKE.STR", DataFrame; delim=' ', ignorerepeated=true, header=["A_H", "Z_H", "TKE", "Y", "σY"], skipto = 2);

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
#####
# Y(A)
function Y_A(dy, A)
    Y = distributie_unidym(Int[], Float64[], Float64[])
    for A_H in minimum(dy.A_H):maximum(dy.A_H)
        # Y(A) = Σ_(Z, TKE) Y(A, Z, TKE)
        # σY(A) = sqrt[Σ_(Z, TKE) σY(A, Z, TKE)^2]
        y_A = sum(dy.Y[dy.A_H .== A_H])
        σ_y_A = sqrt(sum(dy.σY[dy.A_H .== A_H].^2))
        push!(Y.x, A_H)
        push!(Y.y, y_A)
        push!(Y.σ, σ_y_A)
        if A_H != A - A_H
            push!(Y.x, A - A_H)
            push!(Y.y, y_A)
            push!(Y.σ, σ_y_A)
        else
            Y.y[Y.x .== A_H] .+= y_A
            Y.σ[Y.x .== A_H] .+= σ_y_A
        end
    end
    # Normarea distributiei
    f = 200/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f
    return Y
end
# TKE(A)
function TKE_A(dy)
    tke = distributie_unidym(Int[], Float64[], Float64[])
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
        return NaN
    end 
end
# Q(A,Z)
function Q_A_Z(A, Z, df, limInfA_H, limSupA_H)
    Q = distributie_bidym(Int[], Int[], Float64[], Float64[])
    D = df.D[(df.A .== A) .& (df.Z .== Z)][1]
    σ_D = df.σD[(df.A .== A) .& (df.Z .== Z)][1]
    for A_H in limInfA_H:limSupA_H
        for Z_H in minimum(df.Z[df.A .== A_H]):maximum(df.Z[df.A .== A_H])
            A_L = A - A_H
            Z_L = Z - Z_H
            if isassigned(df.D[(df.A .== A_H) .& (df.Z .== Z_H)], 1) && isassigned(df.D[(df.A .== A_L) .& (df.Z .== Z_L)], 1)
                D_H = df.D[(df.A .== A_H) .& (df.Z .== Z_H)][1]
                σ_D_H = df.σD[(df.A .== A_H) .& (df.Z .== Z_H)][1]
                D_L = df.D[(df.A .== A_L) .& (df.Z .== Z_L)][1]
                σ_D_L = df.σD[(df.A .== A_L) .& (df.Z .== Z_L)][1]
                q = D - (D_H + D_L)
                if q > 0
                    # Energiile sunt salvate în MeV
                    push!(Q.y, q *1e-3)
                    push!(Q.σ, sqrt(σ_D^2 + σ_D_H^2 + σ_D_L^2) *1e-3)
                    push!(Q.x_1, A_H)
                    push!(Q.x_2, Z_H)
                end
            end
        end
    end 
    return Q
end
# TXE(A,Z)
function TXE_A_Z(q_A_Z, tke_A, df, A, Z, εₙ)
    txe_A_Z = distributie_bidym(Int[], Int[], Float64[], Float64[])
    Sₙ = Energie_separare(1, 0, A, Z, df)
    for A_H in minimum(q_A_Z.x_1):maximum(q_A_Z.x_1)
        for Z_H in minimum(q_A_Z.x_2[q_A_Z.x_1 .== A_H]):maximum(q_A_Z.x_2[q_A_Z.x_1 .== A_H])
            if isassigned(tke_A.y[tke_A.x .== A_H], 1)
                Q = q_A_Z.y[(q_A_Z.x_1 .== A_H) .& (q_A_Z.x_2 .== Z_H)][1]
                σ_Q = q_A_Z.σ[(q_A_Z.x_1 .== A_H) .& (q_A_Z.x_2 .== Z_H)][1]
                TKE = tke_A.y[tke_A.x .== A_H][1]
                σ_TKE = tke_A.σ[tke_A.x .== A_H][1]
                TXE = Q - TKE + Sₙ[1] + εₙ
                σ_TXE = sqrt(σ_Q^2 + Sₙ[2]^2 + σ_TKE^2)
                if TXE > 0
                    push!(txe_A_Z.x_1, A_H)
                    push!(txe_A_Z.x_2, Z_H)
                    push!(txe_A_Z.y, TXE)
                    push!(txe_A_Z.σ, σ_TXE)
                end
            end
        end
    end
    return txe_A_Z
end
# Calculul T_max
function T_max(A, TXE, C)
    return sqrt(C * TXE/A)
end
# Calculul E_f pentru fragmentele L & H
function E_f(A, A_H, TKE)
    A_L = A - A_H
    return [(A_H/A_L)*TKE/A, (A_L/A_H)*TKE/A]
end
# Calculul u₁ si u₂
function u_1_2(E, E_f, T_max)
    return [(sqrt(E) - sqrt(E_f))^2 /T_max, (sqrt(E) + sqrt(E_f))^2 /T_max]
end
# Calculul integralei E₁(z)
function E_1(z)
    f(x) = exp(-x)/x
    return quadgk(f, z, Inf)[1]
end
# Calculul integralei γ(a,x)
function γ(a, x)
    f(u) = u^(a-1) * exp(-u)
    return quadgk(f, 0, x)[1]
end
# Calculul punctual al lui N(E, E_f)
function N_E_i(E_F, T_MAX, u_1, u_2)
    return (1/3*sqrt(E_F*T_MAX))*(E_1(u_2)*u_2^(3/2) - E_1(u_1)*u_1^(3/2) + γ(3/2, u_2) - γ(3/2, u_1))
end
# Constructia lui N(E)
function N_E(A, Z, E_min, h_E, E_max, limInfA_H, limSupA_H, txe_A_Z, tke_A, y_A)
    n_E = distributie_unidym(Float64[], Float64[], Float64[])
    E = E_min:h_E:E_max
    for i_E in eachindex(E)
        Numarator = 0
        Numitor = 0
        for A_H in limInfA_H:limSupA_H
            if isassigned(y_A.y[y_A.x .== A_H], 1)
                Z_H = round(A_H*Z/A -0.5)
                Y_A = y_A.y[y_A.x .== A_H][1]
                if isassigned(txe_A_Z.y[(txe_A_Z.x_1 .== A_H) .& (txe_A_Z.x_2 .== Z_H)], 1)
                    TXE = txe_A_Z.y[(txe_A_Z.x_1 .== A_H) .& (txe_A_Z.x_2 .== Z_H)][1]
                    TKE = tke_A.y[tke_A.x .== A_H][1]
                    T_MAX = T_max(A, TXE, 10)
                    E_F_L = E_f(A, A_H, TKE)[1]
                    E_F_H = E_f(A, A_H, TKE)[2]
                    u_L = u_1_2(E[i_E], E_F_L, T_MAX)
                    u_H = u_1_2(E[i_E], E_F_H, T_MAX)
                    N_L = N_E_i(E_F_L, T_MAX, u_L[1], u_L[2])
                    N_H = N_E_i(E_F_H, T_MAX, u_H[1], u_H[2])
                    N = 0.5*(N_L + N_H)
                    Numarator += Y_A * N
                    Numitor += Y_A
                end
            end
        end
        if Numitor != 0
            push!(n_E.x, E[i_E])
            push!(n_E.y, Numarator/Numitor)
            #push!(n_E.σ, 0.0)
        end
    end
    return n_E
end
# Normarea spectrului la un spectru Maxwell
function Normare_Maxwell(n_E, T_M)
    Maxwell(E) = (2/sqrt(π)) * T_M^(-3/2) * sqrt(E) * exp(-E/T_M)
    A_calculat = quadgk(Maxwell, first(n_E.x), last(n_E.x))[1]
    A_model = trapz(n_E.x, n_E.y)
    f = A_calculat/A_model
    for i in eachindex(n_E.x)
        n_E.y[i] *= f/Maxwell(n_E.x[i])
    end
    return n_E
end
# Aici se opreste partea de calcul a programului
#####
# Constructia reprezentarilor grafice
function Grafic_plot(distributie, titlu, eticheta, axa_x, axa_y, scalare_inf, scalare_sup, scala_x, scala_y, x_maxim, culoare)
    plt = plot(
        distributie.x, 
        distributie.y, 
        xlims = (minimum(distributie.x), x_maxim),
        ylims = (minimum(distributie.y)*scalare_inf, maximum(distributie.y)*scalare_sup),
        xlabel = axa_x, 
        ylabel = axa_y, 
        framestyle = :box,
        label = eticheta,
        title = titlu,
        xscale = scala_x,
        yscale = scala_y,
        minorgrid = true,
        size = (1000, 1000),
        dpi = 600,
        color = culoare
    )
    return plt
end
function Grafic_afisare(plt, titlu)
    display(plt)
    savefig(plt, "Grafice/T4_$(titlu).png")
end
#####
# Apelarea functiilor definite pentru executia programului
A₀ = 236;
Z₀ = 92;
εₙ = 0;
limInfA_H = 118;
limSupA_H = 160;
E_min = 1e-1;
h_E = E_min/10;
E_max = 20;

y_A = Y_A(dy, A₀);
tke_A = TKE_A(dy);
q_A_Z = Q_A_Z(A₀, Z₀, df, limInfA_H, limSupA_H);
txe_A_Z = TXE_A_Z(q_A_Z, tke_A, df, A₀, Z₀, εₙ);
n_E = N_E(A₀, Z₀, E_min, h_E, E_max, limInfA_H, limSupA_H, txe_A_Z, tke_A, y_A);
# <E_spectru> = 3/2 * T_M_echivalent
T_M = 2/3 * trapz(n_E.x, n_E.y .* n_E.x)/trapz(n_E.x, n_E.y);
n_E = DataFrame(x = n_E.x, y = n_E.y);
n_E_Maxwell = Normare_Maxwell(copy(n_E), T_M);

Plot_n_E_liniar = Grafic_plot(n_E, "Spectrul neutronilor prompți în scală liniară", "", "E [MeV]", "N(E)", 1, 1.1, :identity, :identity, 16, :red);
Grafic_afisare(Plot_n_E_liniar, "Plot_n_E_liniar");

Plot_n_E_logaritmic = Grafic_plot(n_E, "Spectrul neutronilor prompți în scală logaritmică", "", "E [MeV]", "N(E)", 1, 5, :identity, :log10, maximum(n_E.x), :red);
Grafic_afisare(Plot_n_E_logaritmic, "Plot_n_E_logaritmic");

Plot_n_E_Maxwell = Grafic_plot(n_E_Maxwell, "Spectrul neutronilor prompți normat la un spectru Maxwell echivalent", latexstring("\$\\mathrm{T_M} = $(round(T_M, digits = 3))\$ MeV"), "E [MeV]", latexstring("Raportul spectrului neutronilor prompți la distribuția Maxwell având \$\\mathrm{T_M} = $(round(T_M, digits = 3))\$ MeV"), 2, 1.05, :log10, :identity, 10, :red);
hline!(Plot_n_E_Maxwell, [1.0], ls = :dash, label = "", color = :navy);
Grafic_afisare(Plot_n_E_Maxwell, "Plot_n_E_Maxwell");