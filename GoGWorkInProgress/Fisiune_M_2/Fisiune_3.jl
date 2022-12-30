using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul pentru partitionarea energiei totale de excitatie TXE spre fragmentele L & H

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
    for i in minimum(dy.A_H):maximum(dy.A_H)
        # Y(A) = Σ_(Z, TKE) Y(A, Z, TKE)
        # σY(A) = sqrt[Σ_(Z, TKE) σY(A, Z, TKE)^2]
        Suma_Y = sum(dy.Y[dy.A_H .== i])
        Suma_σ = sqrt(sum(dy.σY[dy.A_H .== i].^2))
        push!(Y.x, i)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
        push!(Y.x, A - i)
        push!(Y.y, Suma_Y)
        push!(Y.σ, Suma_σ)
    end
    # Stergerea elementelor duplicate de la A₀/2
    if iseven(A)
        index = findfirst(x -> x == Int(A/2), Y.x)
        deleteat!(Y.x, index)
        deleteat!(Y.y, index)
        deleteat!(Y.σ, index)
    end
    # Normarea distributiei
    f = 200/sum(Y.y)
    Y.y .= Y.y * f
    Y.σ .= Y.σ * f
    return Y
end

function TKE_A(dy)
    TKE = distributie_unidym(Int[],Float64[], Float64[])
    for i in minimum(dy.A_H):maximum(dy.A_H)
        Numarator = 0
        Numitor = sum(dy.Y[dy.A_H .== i])
        Suma_σ² = 0
        # TKE(A) = Σ_(TKE) TKE * Y(A, TKE)/Σ_(TKE) Y(A, TKE)
        # σTKE(A) = [1/Σ_(TKE) Y(A, TKE)] * sqrt[Σ_(TKE) (TKE - TKE(A)) * σY(A, TKE)^2]
        for j in minimum(dy.TKE[dy.A_H .== i]):maximum(dy.TKE[dy.A_H .== i])
            Y_A_TKE = sum(dy.Y[(dy.A_H .== i) .& (dy.TKE .== j)])
            Numarator += j * Y_A_TKE
        end
        tke_A = Numarator/Numitor
        for j in minimum(dy.TKE[dy.A_H .== i]):maximum(dy.TKE[dy.A_H .== i])
            σY_A_TKE = sqrt(sum(dy.σY[(dy.A_H .== i) .& (dy.TKE .== j)].^2))
            Suma_σ² += (j - tke_A)^2 * σY_A_TKE^2
        end
        push!(TKE.x, i)
        push!(TKE.y, tke_A)
        push!(TKE.σ, sqrt(Suma_σ²)/Numitor)
    end  
    return TKE
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

function Sn_A_Z(A, Z, df, limInfA_H, limSupA_H)
    sn = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in limInfA_H:limSupA_H
        if isassigned(df.Z[df.A .== A_H], 1)
            for Z_H in minimum(df.Z[df.A .== A_H]):maximum(df.Z[df.A .== A_H])
                S = Energie_separare(1, 0, A_H, Z_H, df)
                if !isnan(S[1])
                    push!(sn.y, S[1])
                    push!(sn.σ, S[2])
                    push!(sn.x_1, A_H)
                    push!(sn.x_2, Z_H)
                end
                S = Energie_separare(1, 0, A - A_H, Z - Z_H, df)
                if !isnan(S[1])
                    push!(sn.y, S[1])
                    push!(sn.σ, S[2])
                    push!(sn.x_1, A - A_H)
                    push!(sn.x_2, Z - Z_H)
                end
            end
        end
    end
    return sn
end

# Functia de mediere a Sn(A, Z) pe distributia p(A, Z) considerand toate fragmentele
function Sn_A(Sn, A, Z)
    sn = distributie_unidym(Int[], Float64[], Float64[])
    for A_H in minimum(Sn.x_1):maximum(Sn.x_1)
        Numarator = 0
        Numitor = 0
        Sigma_temp² = 0
        Z_UCD = Z*A_H/A
        Z_p = Z_UCD - 0.5
        # Sn(A) = Σ_Z Sn(A, Z)*p(A, Z)/Σ_Z p(A, Z)
        # σSn(A) = sqrt[Σ_Z σSn(A, Z)^2 *p(A, Z)^2]/Σ_Z p(A, Z)
        for j = 1:length(Sn.x_2[Sn.x_1 .== A_H])
            P_A_Z = p_A_Z(Sn.x_2[Sn.x_1 .== A_H][j], Z_p)
            Numarator += Sn.y[Sn.x_1 .== A_H][j] * P_A_Z
            Numitor += P_A_Z
            Sigma_temp² += (P_A_Z * Sn.σ[Sn.x_1 .== A_H][j])^2
        end
        if Numitor != 0
            push!(sn.y, Numarator/Numitor)
            push!(sn.σ, sqrt(Sigma_temp²)/Numitor)
            push!(sn.x, A_H)
        end
    end
    return sn
end

# Distributia izobara de sarcina cu  σ = rms(A) ≈ 0.6 & ΔZₚ ≈ 0.5 (Gaussiana)
function p_A_Z(Z, Z_p)
    return 1/(sqrt(2*π) * 0.6) * exp(-(Z - Z_p)^2 /(2*0.6^2))
end

function Q_A_Z(A, Z, df, limInfA_H, limSupA_H)
    Q = distributie_bidym(Int[], Int[], Float64[], Float64[])
    D = df.D[(df.A .== A) .& (df.Z .== Z)][1]
    σ_D = df.σD[(df.A .== A) .& (df.Z .== Z)][1]
    for A_H in limInfA_H:limSupA_H
        if isassigned(df.Z[df.A .== A_H], 1)
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
    end 
    return Q
end

# Functia de mediere a Q(A, Z) pe distributia p(A, Z) considerand toate fragmentele
function Q_A(Q, A, Z)
    q = distributie_unidym(Int[], Float64[], Float64[])
    for A_H in minimum(Q.x_1):maximum(Q.x_1)
        Numarator = 0
        Numitor = 0
        Sigma_temp² = 0
        Z_UCD = Z*A_H/A
        Z_p = Z_UCD - 0.5
        # Q(A) = Σ_Z Q(A, Z)*p(A, Z)/Σ_Z p(A, Z)
        # σQ(A) = sqrt[Σ_Z σQ(A, Z)^2 *p(A, Z)^2]/Σ_Z p(A, Z)
        for j = 1:length(Q.x_2[Q.x_1 .== A_H])
            P_A_Z = p_A_Z(Q.x_2[Q.x_1 .== A_H][j], Z_p)
            Numarator += Q.y[Q.x_1 .== A_H][j] * P_A_Z
            Numitor += P_A_Z
            Sigma_temp² += (P_A_Z * Q.σ[Q.x_1 .== A_H][j])^2
        end
        if Numitor != 0
            push!(q.y, Numarator/Numitor)
            push!(q.σ, sqrt(Sigma_temp²)/Numitor)
            push!(q.x, A_H)
        end
    end
    return q
end

function TXE_A(q_A, tke_A, df, A, Z, εₙ)
    TXE = distributie_unidym(Int[],Float64[], Float64[])
    Sₙ = Energie_separare(1, 0, A, Z, df)
    for A_H in minimum(q_A.x):maximum(q_A.x)
        if isassigned(tke_A.y[tke_A.x .== A_H], 1) && !isnan(Sₙ[1])
            # TXE(A) = Q(A) - TKE(A) + Sₙ + εₙ
            push!(TXE.x, A_H)
            push!(TXE.y, q_A.y[q_A.x .== A_H][1] + Sₙ[1] + εₙ - tke_A.y[tke_A.x .== A_H][1])
            push!(TXE.σ, sqrt(q_A.σ[q_A.x .== A_H][1]^2 + Sₙ[2]^2 + tke_A.σ[tke_A.x .== A_H][1]^2))
        end
    end  
    return TXE
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

# Apelarea functiilor definite pentru executia programului
A₀ = 236
Z₀ = 92
εₙ = 0
limInfA_H = 118
limSupA_H = 170

y_A = Sortare_distributie(Y_A(dy, A₀))
tke_A = Sortare_distributie(TKE_A(dy))
q_A = Sortare_distributie(Q_A(Q_A_Z(A₀, Z₀, df, limInfA_H, limSupA_H), A₀, Z₀))
sn_A = Sortare_distributie(Sn_A(Sn_A_Z(A₀, Z₀, df, limInfA_H, limSupA_H), A₀, Z₀))
txe_A = Sortare_distributie(TXE_A(q_A, tke_A, df, A₀, Z₀, εₙ))

Q_med = Medie_distributie(q_A, y_A, firstindex(q_A.x), lastindex(q_A.x))
mid_index_sn_A = Int((length(sn_A.x) + 1 )/2);
Sn_H_med = Medie_distributie(sn_A, y_A, firstindex(sn_A.x), mid_index_sn_A)
Sn_L_med = Medie_distributie(sn_A, y_A, mid_index_sn_A, lastindex(sn_A.x))
TXE_med = Medie_distributie(txe_A, y_A, firstindex(txe_A.x), lastindex(txe_A.x))