using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul pentru partitionarea energiei totale de excitatie TXE intre fragmentele L & H folosind modelarea la sciziune

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

# Citire fisiere de date
dm = CSV.read("Data_files/Defecte_masa/AUDI2021.csv", DataFrame; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σD"]);
dy = CSV.read("Data_files/Yield/U5YAZTKE.STR", DataFrame; delim=' ', ignorerepeated=true, header=["A_H", "Z_H", "TKE", "Y", "σY"], skipto = 2);
dβ₀ = CSV.read("Data_files/Parametrizari_auxiliare/B2MOLLER.ANA", DataFrame; delim=' ', ignorerepeated=true, header=["Z", "A", "β"], skipto = 2);
dGC = CSV.read("Data_files/Parametrizari_auxiliare/SZSN.GC", DataFrame; delim=' ', ignorerepeated=true, header=["n", "S_N", "S_Z"], skipto = 2);
dν_Gook = CSV.read("Data_files/Date_experimentale/Multiplicitate_n/U5NUAGOOK.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["A", "ν", "σν"], skipto = 2);
dν_Maslin = CSV.read("Data_files/Date_experimentale/Multiplicitate_n/U5NUAMASLIN.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["A", "ν", "σν"], skipto = 2);
dν_Nishio = CSV.read("Data_files/Date_experimentale/Multiplicitate_n/U5NUANISHIO.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["A", "ν", "σν"], skipto = 2);
dν_Vorobyev = CSV.read("Data_files/Date_experimentale/Multiplicitate_n/U5NUAVORO.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["A", "ν", "σν"], skipto = 2);

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
# Calculul energiei de separare a particulei (A_part,Z_part) din nucleul (A,Z)
function Energie_separare(A_part, Z_part, A, Z, dm)
    # Verificarea existentei radionuclizilor folositi in libraria de date
    if isassigned(dm.D[(dm.A .== A_part) .& (dm.Z .== Z_part)], 1) && isassigned(dm.D[(dm.A .== A) .& (dm.Z .== Z)], 1) && isassigned(dm.D[(dm.A .== A - A_part) .& (dm.Z .== Z - Z_part)], 1)
        D_part = dm.D[(dm.A .== A_part) .& (dm.Z .== Z_part)][1]
        σ_part = dm.σD[(dm.A .== A_part) .& (dm.Z .== Z_part)][1]
        D = dm.D[(dm.A .== A - A_part) .& (dm.Z .== Z - Z_part)][1]
        σᴰ = dm.σD[(dm.A .== A - A_part) .& (dm.Z .== Z - Z_part)][1]    
        S = (D + D_part - dm.D[(dm.A .== A) .& (dm.Z .== Z)][1])*1e-3
        σˢ = sqrt(σ_part^2 + σᴰ^2 + (dm.σD[(dm.A .== A) .& (dm.Z .== Z)][1])^2)*1e-3
        return [S, σˢ]
    else 
        return NaN
    end 
end
# Calculul Sₙ(A,Z) pentru fragmentele L & H
function Sn_A_Z(A, Z, dm, limInfA_H, limSupA_H)
    sn = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in limInfA_H:limSupA_H
        for Z_H in minimum(dm.Z[dm.A .== A_H]):maximum(dm.Z[dm.A .== A_H])
            S_H = Energie_separare(1, 0, A_H, Z_H, dm)
            S_L = Energie_separare(1, 0, A - A_H, Z - Z_H, dm)
            if !isnan(S_H[1]) && !isnan(S_L[1])
                push!(sn.y, S_H[1])
                push!(sn.σ, S_H[2])
                push!(sn.x_1, A_H)
                push!(sn.x_2, Z_H)
                if A - A_H != A_H
                    push!(sn.y, S_L[1])
                    push!(sn.σ, S_L[2])
                    push!(sn.x_1, A - A_H)
                    push!(sn.x_2, Z - Z_H)
                end
            end
        end
    end
    return sn
end
# Distributia izobara de sarcina cu σ = rms(A) ≈ 0.6 & ΔZₚ ≈ 0.5 (Gaussiana)
function p_A_Z(Z, Z_p)
    return 1/(sqrt(2*π) * 0.6) * exp(-(Z - Z_p)^2 /(2*0.6^2))
end
# Q(A,Z)
function Q_A_Z(A, Z, dm, limInfA_H, limSupA_H)
    Q = distributie_bidym(Int[], Int[], Float64[], Float64[])
    D = dm.D[(dm.A .== A) .& (dm.Z .== Z)][1]
    σ_D = dm.σD[(dm.A .== A) .& (dm.Z .== Z)][1]
    for A_H in limInfA_H:limSupA_H
        for Z_H in minimum(dm.Z[dm.A .== A_H]):maximum(dm.Z[dm.A .== A_H])
            A_L = A - A_H
            Z_L = Z - Z_H
            if isassigned(dm.D[(dm.A .== A_H) .& (dm.Z .== Z_H)], 1) && isassigned(dm.D[(dm.A .== A_L) .& (dm.Z .== Z_L)], 1)
                D_H = dm.D[(dm.A .== A_H) .& (dm.Z .== Z_H)][1]
                σ_D_H = dm.σD[(dm.A .== A_H) .& (dm.Z .== Z_H)][1]
                D_L = dm.D[(dm.A .== A_L) .& (dm.Z .== Z_L)][1]
                σ_D_L = dm.σD[(dm.A .== A_L) .& (dm.Z .== Z_L)][1]
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
function TXE_A_Z(q_A_Z, tke_A, dm, A, Z, εₙ)
    txe_A_Z = distributie_bidym(Int[], Int[], Float64[], Float64[])
    Sₙ = Energie_separare(1, 0, A, Z, dm)
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
function Functie_liniara(x, a, b)
    return a*x + b
end
# Parametrizarea factorului de deformare la sciziune
function Beta_sciziune()
    β_sciz = distributie_unidym(Int[], Float64[], Float64[])
    a = 0.58/13
    b = -28*a
    for Z in 28:41
        push!(β_sciz.x, Z)
        push!(β_sciz.y, Functie_liniara(Z, a, b))
        push!(β_sciz.σ, 0.0)
    end
    for Z in 42:44
        push!(β_sciz.x, Z)
        push!(β_sciz.y, 0.58)
        push!(β_sciz.σ, 0.0)
    end
    a = -0.58/6
    b = -50*a
    for Z in 45:50
        push!(β_sciz.x, Z)
        push!(β_sciz.y, Functie_liniara(Z, a, b))
        push!(β_sciz.σ, 0.0)
    end
    a = 0.6/15
    b = -50*a
    for Z in 51:65
        push!(β_sciz.x, Z)
        push!(β_sciz.y, Functie_liniara(Z, a, b))
        push!(β_sciz.σ, 0.0)
    end
    return β_sciz
end
# Calculul punctual al energiei nucleului cu modelul picatura de lichid
function E_LDM(β, A, Z)
    η = (A - 2*Z)/A
    χ = 1 - 1.7826 * η^2
    α² = (5 * β^2)/(4*π)
    X_νs = -χ*(15.4941*A - 17.9439*A^(2/3) * (1 + 0.4*α²))
    X_e = Z^2 * (0.7053 * (1 - 0.2*α²)/(A^(1/3)) - 1.1529/A)
    return X_νs + X_e
end
# Calculul energiei de extra-deformare pentru toate fragmentele
function ΔE_deformare_A_Z(A, Z, dβ₀, β_sciz, limInfA_H, limSupA_H)
    ΔE = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in limInfA_H:limSupA_H
        # Verificarea existentei unei valori pentru β_gs in baza de date Moller & Nix
        if isassigned(dβ₀.Z[dβ₀.A .== A_H], 1)
            for Z_H in minimum(dβ₀.Z[dβ₀.A .== A_H]):maximum(dβ₀.Z[dβ₀.A .== A_H])
                A_L = A - A_H
                Z_L = Z - Z_H
                #Verificarea existentei valorilor β_sciziune din parametrizare
                if isassigned(β_sciz.y[β_sciz.x .== Z_H], 1) && isassigned(β_sciz.y[β_sciz.x .== Z_L], 1)
                    if isassigned(dβ₀.β[(dβ₀.A .== A_L) .& (dβ₀.Z .== Z_L)], 1)
                        E_LDM_H_0 = E_LDM(0, A_H, Z_H)
                        E_LDM_L_0 = E_LDM(0, A_L, Z_L)
                        E_def_sciz_H = E_LDM(β_sciz.y[β_sciz.x .== Z_H][1], A_H, Z_H) - E_LDM_H_0
                        E_def_sciz_L = E_LDM(β_sciz.y[β_sciz.x .== Z_L][1], A_L, Z_L) - E_LDM_L_0
                        E_def_at_H = E_LDM(dβ₀.β[(dβ₀.A .== A_H) .& (dβ₀.Z .== Z_H)][1], A_H, Z_H) - E_LDM_H_0
                        E_def_at_L = E_LDM(dβ₀.β[(dβ₀.A .== A_L) .& (dβ₀.Z .== Z_L)][1], A_L, Z_L) - E_LDM_L_0
                        push!(ΔE.x_1, A_H)
                        push!(ΔE.x_2, Z_H)
                        push!(ΔE.y, abs(E_def_sciz_H - E_def_at_H))
                        push!(ΔE.σ, 0.0)
                        if A_L != A_H
                            push!(ΔE.x_1, A_L)
                            push!(ΔE.x_2, Z_L)
                            push!(ΔE.y, abs(E_def_sciz_L - E_def_at_L))
                            push!(ΔE.σ, 0.0)
                        end
                    end
                end
            end
        end
    end
    return ΔE
end
# Calculul parametrului densitatilor de nivele, sistematica Gilbert-Cameron
function a_Gilbert_Cameron(dGC, A, Z, limInfA_H, limSupA_H)
    a = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in limInfA_H:limSupA_H
        for Z_H in minimum(dGC.n):maximum(dGC.n)
            N_H = A_H - Z_H
            A_L = A - A_H
            Z_L = Z - Z_H
            N_L = A_L - Z_L
            if isassigned(dGC.S_Z[dGC.n .== Z_H], 1) && isassigned(dGC.S_N[dGC.n .== N_H], 1)
                if isassigned(dGC.S_Z[dGC.n .== Z_L], 1) && isassigned(dGC.S_N[dGC.n .== N_L], 1)
                    a_H = A_H * (0.00917*(dGC.S_Z[dGC.n .== Z_H][1] + dGC.S_N[dGC.n .== N_H][1]) + 0.142)
                    a_L = A_L * (0.00917*(dGC.S_Z[dGC.n .== Z_L][1] + dGC.S_N[dGC.n .== N_L][1]) + 0.142)
                    if a_H > 0 && a_L > 0
                        push!(a.x_1, A_H)
                        push!(a.x_2, Z_H)
                        push!(a.y, a_H)
                        push!(a.σ, 0.0)
                        if A_L != A_H
                            push!(a.x_1, A_L)
                            push!(a.x_2, Z_L)
                            push!(a.y, a_L)
                            push!(a.σ, 0.0)
                        end
                    end
                end
            end
        end
    end
    return a
end
# Calculul si partitionarea energiei de excitatie disponibile la sciziune
function E_sciziune(txe, ΔE_def, a, A, Z, limInfA_H, limSupA_H)
    E_sciz = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in limInfA_H:limSupA_H
        A_L = A - A_H
        # Baleiajul fragmentarilor H se realizeaza dupa valorile TXE(A, Z)
        if isassigned(txe.x_2[txe.x_1 .== A_H], 1)
            for Z_H in minimum(txe.x_2[txe.x_1 .== A_H]):maximum(txe.x_2[txe.x_1 .== A_H])
                Z_L = Z - Z_H
                # Verificam existenta valorilor pentru ΔE_deformare ale ambelor fragmente din perechea L-H
                if isassigned(ΔE_def.y[(ΔE_def.x_1 .== A_H) .& (ΔE_def.x_2 .== Z_H)], 1) && isassigned(ΔE_def.y[(ΔE_def.x_1 .== A_L) .& (ΔE_def.x_2 .== Z_L)], 1)
                    # Analog pentru existenta parametrilor densitatilor de nivele
                    if isassigned(a.y[(a.x_1 .== A_H) .& (a.x_2 .== Z_H)], 1) && isassigned(a.y[(a.x_1 .== A_L) .& (a.x_2 .== Z_L)], 1)
                        TXE = txe.y[(txe.x_1 .== A_H) .& (txe.x_2 .== Z_H)][1]
                        ΔE_def_L = ΔE_def.y[(ΔE_def.x_1 .== A_L) .& (ΔE_def.x_2 .== Z_L)][1]
                        ΔE_def_H = ΔE_def.y[(ΔE_def.x_1 .== A_H) .& (ΔE_def.x_2 .== Z_H)][1]
                        ε_sciz = TXE - (ΔE_def_L + ΔE_def_H)
                        a_L = a.y[(a.x_1 .== A_L) .& (a.x_2 .== Z_L)][1]
                        a_H = a.y[(a.x_1 .== A_H) .& (a.x_2 .== Z_H)][1]
                        r = a_L/a_H
                        if ε_sciz > 0
                            push!(E_sciz.x_1, A_H)
                            push!(E_sciz.x_2, Z_H)
                            push!(E_sciz.y, ε_sciz/(1 + r))
                            push!(E_sciz.σ, 0.0)
                            if A_L != A_H
                                push!(E_sciz.x_1, A_L)
                                push!(E_sciz.x_2, Z_L)
                                push!(E_sciz.y, ε_sciz * r/(1 + r))
                                push!(E_sciz.σ, 0.0)
                            end     
                        end
                    end
                end
            end
        end
    end  
    return E_sciz
end
# Calculul energiei de excitatie a fragmentelor total accelerate
function E_excitatie_at(ΔE_def, E_sciz, A, Z, limInfA_H, limSupA_H)
    E_excitatie = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in limInfA_H:limSupA_H
        A_L = A - A_H
        if isassigned(E_sciz.x_2[E_sciz.x_1 .== A_H], 1)
            for Z_H in minimum(E_sciz.x_2[E_sciz.x_1 .== A_H]):maximum(E_sciz.x_2[E_sciz.x_1 .== A_H])
                Z_L = Z - Z_H
                ΔE_def_L = ΔE_def.y[(ΔE_def.x_1 .== A_L) .& (ΔE_def.x_2 .== Z_L)][1]
                ΔE_def_H = ΔE_def.y[(ΔE_def.x_1 .== A_H) .& (ΔE_def.x_2 .== Z_H)][1]
                E_sciz_L = E_sciz.y[(E_sciz.x_1 .== A_L) .& (E_sciz.x_2 .== Z_L)][1]
                E_sciz_H = E_sciz.y[(E_sciz.x_1 .== A_H) .& (E_sciz.x_2 .== Z_H)][1]
                E_excitatie_L = ΔE_def_L + E_sciz_L
                E_excitatie_H = ΔE_def_H + E_sciz_H
                push!(E_excitatie.x_1, A_H)
                push!(E_excitatie.x_2, Z_H)
                push!(E_excitatie.y, E_excitatie_H)
                push!(E_excitatie.σ, 0.0)
                if A_L != A_H
                    push!(E_excitatie.x_1, A_L)
                    push!(E_excitatie.x_2, Z_L)
                    push!(E_excitatie.y, E_excitatie_L)
                    push!(E_excitatie.σ, 0.0)
                end
            end
        end
    end
    return E_excitatie
end
# Calculul raportului final de partitionare R
function R_partitionare(E_excitatie, txe)
    R = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in minimum(txe.x_1):maximum(txe.x_1)
        for Z_H in minimum(txe.x_2[txe.x_1 .== A_H]):maximum(txe.x_2[txe.x_1 .== A_H])
            if isassigned(E_excitatie.y[(E_excitatie.x_1 .== A_H) .& (E_excitatie.x_2 .== Z_H)], 1)
                TXE = txe.y[(txe.x_1 .== A_H) .& (txe.x_2 .== Z_H)][1]
                E_excitatie_H = E_excitatie.y[(E_excitatie.x_1 .== A_H) .& (E_excitatie.x_2 .== Z_H)][1]
                push!(R.x_1, A_H)
                push!(R.x_2, Z_H)
                push!(R.y, E_excitatie_H/TXE)
                push!(R.σ, 0.0)
            end
        end
    end
    return R
end
# Calculul multiplicitatii neutronice prompte a perechilor de fragmente
function ν_pereche(txe, sn, a, A, Z)
    ν = distributie_bidym(Int[], Int[], Float64[], Float64[])
    p = 6.71 - Z^2 * 0.156/A
    q = 0.75 + Z^2 * 0.088/A
    for A_H in minimum(txe.x_1):maximum(txe.x_1)
        A_L = A - A_H
        for Z_H in minimum(txe.x_2[txe.x_1 .== A_H]):maximum(txe.x_2[txe.x_1 .== A_H])
            Z_L = Z - Z_H
            a_L = a.y[(a.x_1 .== A_L) .& (a.x_2 .== Z_L)][1]
            a_H = a.y[(a.x_1 .== A_H) .& (a.x_2 .== Z_H)][1]
            a_med = a_L + a_H
            TXE = txe.y[(txe.x_1 .== A_H) .& (txe.x_2 .== Z_H)][1]
            Tₘ = sqrt(TXE/a_med)
            ε_med = 4*Tₘ/3
            sn_L = sn.y[(sn.x_1 .== A_L) .& (sn.x_2 .== Z_L)][1]
            sn_H = sn.y[(sn.x_1 .== A_H) .& (sn.x_2 .== Z_H)][1]
            Sₙ_med = 0.5*(sn_L + sn_H)
            if TXE >= Sₙ_med
                push!(ν.x_1, A_H)
                push!(ν.x_2, Z_H)
                push!(ν.y, (TXE - q)/(ε_med + Sₙ_med + p))
                push!(ν.σ, 0.0)
            end
        end
    end
    return ν
end
# Calculul multiplicitatii neutronice prompte partitionata pe fiecare fragmentarilor
function ν_partitionat(ν_pereche, R, A, Z)
    ν = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in minimum(R.x_1):maximum(R.x_1)
        A_L = A - A_H
        for Z_H in minimum(R.x_2[R.x_1 .== A_H]):maximum(R.x_2[R.x_1 .== A_H])
            if isassigned(ν_pereche.y[(ν_pereche.x_1 .== A_H) .& (ν_pereche.x_2 .== Z_H)], 1)
                Z_L = Z - Z_H
                mu_pereche = ν_pereche.y[(ν_pereche.x_1 .== A_H) .& (ν_pereche.x_2 .== Z_H)][1]
                R_pereche = R.y[(R.x_1 .== A_H) .& (R.x_2 .== Z_H)][1]
                push!(ν.x_1, A_H)
                push!(ν.x_2, Z_H)
                push!(ν.y, R_pereche * mu_pereche)
                push!(ν.σ, 0.0)
                if A_L != A_H
                    push!(ν.x_1, A_L)
                    push!(ν.x_2, Z_L)
                    push!(ν.y, (1 - R_pereche) * mu_pereche)
                    push!(ν.σ, 0.0)
                end
            end
        end
    end
    return ν
end
# Indicele de mijloc al abscisei unei distributii unidym returnat ca Int
function Indice_mijloc(distributie)
    L = length(distributie.x)
    if iseven(L)
        return Int(L/2)
    else
        return Int((L+1)/2)
    end
end
# Sortarea crescatoare dupa ordonata a structurilor de tip obiect unidimensional
# Si stergerea elementelor duplicate dupa abscisa
function Sortare_distributie(distributie)
    #index_true = unique(i -> distributie.x[i], eachindex(distributie.x))
    #index_delete = setdiff(eachindex(distributie.x), index_true)
    #deleteat!(distributie.x, index_delete)
    #deleteat!(distributie.y, index_delete)
    #deleteat!(distributie.σ, index_delete)
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
# Valoare medie & incertitudine pentru distributii mediate folosind o distributie de yield
function Medie_distributie_yield(distributie, Y, index_min, index_max)
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
                Suma_σ² += (distributie.y[i] - Media_distributiei)^2 * Y.σ[Y.x .== distributie.x[i]][1]^2
            end
        end    
        return [round(Media_distributiei, digits = 3), round(sqrt(Suma_σ²)/Numitor, digits = 5)]    
    else 
        return NaN
    end
end
# Functia de mediere a unei distributii_bidym(A,Z) pe distributia p(A,Z)
function Mediere_distributie_p(distributie_A_Z, A, Z, limInfA_H, limSupA_H)
    distributie_A = distributie_unidym(Int[], Float64[], Float64[])
    # Baleiaj interval fragmente H
    for A_H in limInfA_H:limSupA_H
        Numarator = 0
        Numitor = 0
        Sigma_temp² = 0
        Z_UCD = Z*A_H/A
        Z_p = Z_UCD - 0.5
        for j = 1:length(distributie_A_Z.x_2[distributie_A_Z.x_1 .== A_H])
            P_A_Z = p_A_Z(distributie_A_Z.x_2[distributie_A_Z.x_1 .== A_H][j], Z_p)
            Numarator += distributie_A_Z.y[distributie_A_Z.x_1 .== A_H][j] * P_A_Z
            Numitor += P_A_Z
            Sigma_temp² += (P_A_Z * distributie_A_Z.σ[distributie_A_Z.x_1 .== A_H][j])^2
        end
        if Numitor != 0
            push!(distributie_A.y, Numarator/Numitor)
            push!(distributie_A.σ, sqrt(Sigma_temp²)/Numitor)
            push!(distributie_A.x, A_H)
        end
    end
    # Baleiaj interval fragmente L
    limInfA_L = A - limSupA_H
    limSupA_L = A - limInfA_H
    if limSupA_L == limInfA_H
        limSupA_L -= 1
    end
    for A_L in limInfA_L:limSupA_L
        Numarator = 0
        Numitor = 0
        Sigma_temp² = 0
        Z_UCD = Z*A_L/A
        Z_p = Z_UCD + 0.5
        for j = 1:length(distributie_A_Z.x_2[distributie_A_Z.x_1 .== A_L])
            P_A_Z = p_A_Z(distributie_A_Z.x_2[distributie_A_Z.x_1 .== A_L][j], Z_p)
            Numarator += distributie_A_Z.y[distributie_A_Z.x_1 .== A_L][j] * P_A_Z
            Numitor += P_A_Z
            Sigma_temp² += (P_A_Z * distributie_A_Z.σ[distributie_A_Z.x_1 .== A_L][j])^2
        end
        if Numitor != 0
            push!(distributie_A.y, Numarator/Numitor)
            push!(distributie_A.σ, sqrt(Sigma_temp²)/Numitor)
            push!(distributie_A.x, A_L)
        end
    end
    return distributie_A
end
# Trecerea datelor experimentale de la tot intervalul A la intervalul fragmentelor grele
function ν_pair_experimental(dν, A, limInfA_H, limSupA_H)
    ν_pair_exp = distributie_unidym(Int[], Float64[], Float64[])
    for A_H in limInfA_H:limSupA_H
        A_L = A - A_H
        if isassigned(dν.A[dν.A .== A_H], 1) && isassigned(dν.A[dν.A .== A_L], 1)
            push!(ν_pair_exp.x, A_H)
            push!(ν_pair_exp.y, dν.ν[dν.A .== A_H][1] + dν.ν[dν.A .== A_L][1])
            push!(ν_pair_exp.σ, sqrt(dν.σν[dν.A .== A_H][1]^2 + dν.σν[dν.A .== A_L][1]^2))
        end
    end
    return ν_pair_exp
end
# Calculul raportului νₕ/ν_pair
function Raport_νH_νPair(ν, ν_pereche, A)
    raport = distributie_unidym(Int[], Float64[], Float64[])
    for i in eachindex(ν_pereche.x)
        A_H = ν_pereche.x[i]
        if isassigned(ν.x[ν.x .== A_H], 1)
            ν_p = ν_pereche.y[i]
            ν_H = ν.y[ν.x .== A_H][1]
            f = ν_H/ν_p
            push!(raport.x, A_H)
            push!(raport.y, f)
            σν_p = ν_pereche.σ[i]
            if σν_p != 0
                push!(raport.σ, (1/ν_p) * sqrt(ν.σ[ν.x .== A_H][1]^2 + f^2 * σν_p^2))
            else
                push!(raport.σ, 0.0)
            end
        end
    end
    return raport
end
# Aici se opreste partea de calcul a programului
#####
# Constructia reprezentarilor grafice
function Grafic_scatter(x, y, σ, titlu, eticheta, culoare, forma, marime, axa_x, axa_y, scalare_inf, scalare_sup, x_minim, x_maxim)
    plt = scatter(
        x, 
        y,
        yerr = σ, 
        ylims = (minimum(y)*scalare_inf, maximum(y)*scalare_sup),
        xlims = (x_minim, x_maxim),
        xlabel = axa_x, 
        ylabel = axa_y, 
        framestyle = :box,
        label = eticheta,
        title = titlu,
        minorgrid = true,
        size = (1000, 950),
        dpi = 600,
        color = culoare,
        shape = forma,
        markersize = marime
    )
    return plt
end
function Grafic_scatter(x, y, σ, eticheta, culoare, forma, marime, plt)
    scatter!(plt,
        x, 
        y,
        yerr = σ, 
        label = eticheta,
        color = culoare,
        shape = forma,
        markersize = marime
    )
    return plt
end
function Grafic_unire_linie(x, y, σ, plt, eticheta, culoare)
    plot!(plt, 
    x, 
    y,
    ribbon = σ,
    #fillalpha = .3,
    label = eticheta,
    color = culoare
    )
    return plt
end
function Grafic_textbox(x, y, plt, text)
    annotate!(plt, 
    x, 
    y, 
    text
    )
    return plt
end
function Grafic_linie_medie_vertical(plt, valoare_medie, eticheta, culoare)
    vline!(plt, 
    [valoare_medie], 
    ls = :dashdot, 
    label = eticheta,
    color = culoare
    )
    return plt
end
function Grafic_linie_medie_orizontal(plt, valoare_medie, eticheta, culoare)
    hline!(plt, 
    [valoare_medie], 
    ls = :dashdot, 
    label = eticheta,
    color = culoare
    )
    return plt
end
function Grafic_afisare(plt, titlu)
    display(plt)
    savefig(plt, "Grafice/T3_$(titlu).png")
end
#####
# Apelarea functiilor definite pentru executia programului
A₀ = 236;
Z₀ = 92;
εₙ = 0;
limInfA_H = 118;
limSupA_H = 160;

y_A = Sortare_distributie(Y_A(dy, A₀));
tke_A = Sortare_distributie(TKE_A(dy));
q_A_Z = Q_A_Z(A₀, Z₀, dm, limInfA_H, limSupA_H);
q_A = Sortare_distributie(Mediere_distributie_p(q_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
sn_A_Z = Sn_A_Z(A₀, Z₀, dm, limInfA_H, limSupA_H);
sn_A = Sortare_distributie(Mediere_distributie_p(sn_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
txe_A_Z = TXE_A_Z(q_A_Z, tke_A, dm, A₀, Z₀, εₙ);
txe_A = Sortare_distributie(Mediere_distributie_p(txe_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
β_sciz = Beta_sciziune();
ΔE_def_A_Z = ΔE_deformare_A_Z(A₀, Z₀, dβ₀, β_sciz, limInfA_H, limSupA_H);
ΔE_def_A = Sortare_distributie(Mediere_distributie_p(ΔE_def_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
a_A_Z = a_Gilbert_Cameron(dGC, A₀, Z₀, limInfA_H, limSupA_H);
a_A = Sortare_distributie(Mediere_distributie_p(a_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
E_sciz_A_Z = E_sciziune(txe_A_Z, ΔE_def_A_Z, a_A_Z, A₀, Z₀, limInfA_H, limSupA_H);
E_sciz_A = Sortare_distributie(Mediere_distributie_p(E_sciz_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
E_excitatie_A_Z = E_excitatie_at(ΔE_def_A_Z, E_sciz_A_Z, A₀, Z₀, limInfA_H, limSupA_H);
E_excitatie_A = Sortare_distributie(Mediere_distributie_p(E_excitatie_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
R_A_Z = R_partitionare(E_excitatie_A_Z, txe_A_Z);
R_A = Sortare_distributie(Mediere_distributie_p(R_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
ν_pereche_A_Z = ν_pereche(txe_A_Z, sn_A_Z, a_A_Z, A₀, Z₀);
ν_pereche_A = Sortare_distributie(Mediere_distributie_p(ν_pereche_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
ν_A_Z = ν_partitionat(ν_pereche_A_Z, R_A_Z, A₀, Z₀);
ν_A = Sortare_distributie(Mediere_distributie_p(ν_A_Z, A₀, Z₀, limInfA_H, limSupA_H));
raport_νH_νPair_Model = Raport_νH_νPair(ν_A, ν_pereche_A, A₀);

ν_pereche_Gook = ν_pair_experimental(copy(dν_Gook), A₀, limInfA_H, limSupA_H);
ν_Gook = DataFrame(x = dν_Gook.A, y = dν_Gook.ν, σ = dν_Gook.σν);
raport_νH_νPair_Gook = Raport_νH_νPair(ν_Gook, ν_pereche_Gook, A₀);

ν_pereche_Maslin = ν_pair_experimental(copy(dν_Maslin), A₀, limInfA_H, limSupA_H);
ν_Maslin = DataFrame(x = dν_Maslin.A, y = dν_Maslin.ν, σ = dν_Maslin.σν);
raport_νH_νPair_Maslin = Raport_νH_νPair(ν_Maslin, ν_pereche_Maslin, A₀);

ν_pereche_Nishio = ν_pair_experimental(copy(dν_Nishio), A₀, limInfA_H, limSupA_H);
ν_Nishio = DataFrame(x = dν_Nishio.A, y = dν_Nishio.ν, σ = dν_Nishio.σν);
raport_νH_νPair_Nishio = Raport_νH_νPair(ν_Nishio, ν_pereche_Nishio, A₀);

ν_pereche_Vorobyev = ν_pair_experimental(copy(dν_Vorobyev), A₀, limInfA_H, limSupA_H);
ν_Vorobyev = DataFrame(x = dν_Vorobyev.A, y = dν_Vorobyev.ν, σ = dν_Vorobyev.σν);
raport_νH_νPair_Vorobyev = Raport_νH_νPair(ν_Vorobyev, ν_pereche_Vorobyev, A₀);


Q_med = Medie_distributie_yield(q_A, y_A, firstindex(q_A.x), lastindex(q_A.x));
Plot_Q_A = Grafic_scatter(q_A.x, q_A.y, q_A.σ, "Q(A) obținut prin medierea Q(A,Z) al fragmentelor pe distribuția izobară de sarcină", "", :orangered, :circle, 6, L"\mathrm{A_H}", "Q [MeV]", 0.98, 1.03, first(q_A.x)-0.5, last(q_A.x)+0.5);
Plot_Q_A = Grafic_textbox(151.5, 200, Plot_Q_A, latexstring("\$<\\mathrm{Q}> =\$ $(round(Q_med[1], digits = 3)) \$\\pm\$ $(round(Q_med[2], digits = 5)) MeV"));
Plot_Q_A = Grafic_unire_linie(q_A.x, q_A.y, q_A.σ, Plot_Q_A, "", :orangered);
Grafic_afisare(Plot_Q_A, "Q_A");

Sn_H_med = Medie_distributie_yield(sn_A, y_A, firstindex(sn_A.x), Indice_mijloc(sn_A));
Sn_L_med = Medie_distributie_yield(sn_A, y_A, Indice_mijloc(sn_A), lastindex(sn_A.x));
Plot_Sn_A = Grafic_scatter(sn_A.x, sn_A.y, sn_A.σ, "Energia de separare a neutronului din fragmentele de fisiune", "", :orangered, :circle, 6, "A", latexstring("\$\\mathrm{S_n}\$ [MeV]"), 0.95, 1.05, first(sn_A.x)-0.5, last(sn_A.x)+0.5);
Plot_Sn_A = Grafic_textbox(100, 7.25, Plot_Sn_A, latexstring("\$<\\mathrm{S_n^H}> =\$ $(round(Sn_H_med[1], digits = 3)) \$\\pm\$ $(round(Sn_H_med[2], digits = 5)) MeV"));
Plot_Sn_A = Grafic_textbox(100, 7.1, Plot_Sn_A, latexstring("\$<\\mathrm{S_n^L}> =\$ $(round(Sn_L_med[1], digits = 3)) \$\\pm\$ $(round(Sn_L_med[2], digits = 5)) MeV"));
Plot_Sn_A = Grafic_unire_linie(sn_A.x, sn_A.y, sn_A.σ, Plot_Sn_A, "", :orangered);
Grafic_afisare(Plot_Sn_A, "Sn_A");

TXE_med = Medie_distributie_yield(txe_A, y_A, firstindex(txe_A.x), lastindex(txe_A.x));
Plot_TXE_A = Grafic_scatter(txe_A.x, txe_A.y, txe_A.σ, "TXE(A) obținut prin medierea TXE(A,Z) al fragmentelor pe distribuția izobară de sarcină", "", :orangered, :circle, 5, L"\mathrm{A_H}", "TXE [MeV]", 0.95, 1.05, first(txe_A.x)-0.5, last(txe_A.x)+0.5);
Plot_TXE_A = Grafic_textbox(140, 41.5, Plot_TXE_A, latexstring("\$<\\mathrm{TXE}> =\$ $(round(TXE_med[1], digits = 3)) \$\\pm\$ $(round(TXE_med[2], digits = 5)) MeV"));
Plot_TXE_A = Grafic_linie_medie_orizontal(Plot_TXE_A, TXE_med[1], L"<\mathrm{TXE}>", :olive);
Plot_TXE_A = Grafic_unire_linie(txe_A.x, txe_A.y, txe_A.σ, Plot_TXE_A, "", :orangered);
Grafic_afisare(Plot_TXE_A, "TXE_A");

a_H_med = Medie_distributie_yield(a_A, y_A, firstindex(a_A.x), Indice_mijloc(a_A));
a_L_med = Medie_distributie_yield(a_A, y_A, Indice_mijloc(a_A), lastindex(a_A.x));
Plot_a_A = Grafic_scatter(a_A.x, a_A.y, a_A.σ, "Parametrul densității de nivele energetice, sistematica Gilbert-Cameron", "", :orangered, :circle, 6, "A", L"\mathrm{a \;\; [MeV^{-1}]}", 0.95, 1.05, first(a_A.x)-0.5, last(a_A.x)+0.5);
Plot_a_A = Grafic_textbox(a_A.x[Indice_mijloc(a_A)]*0.9, maximum(a_A.y)*1.02, Plot_a_A, latexstring("\$<\\mathrm{a^H}> =\$ $(round(a_H_med[1], digits = 3)) \$\\pm\$ $(round(a_H_med[2], digits = 5)) MeV\$^{-1}\$"));
Plot_a_A = Grafic_textbox(a_A.x[Indice_mijloc(a_A)]*0.9, maximum(a_A.y)*0.99, Plot_a_A, latexstring("\$<\\mathrm{a^L}> =\$ $(round(a_L_med[1], digits = 3)) \$\\pm\$ $(round(a_L_med[2], digits = 5)) MeV\$^{-1}\$"));
Grafic_afisare(Plot_a_A, "a_A");

Plot_E_A = Grafic_scatter(ΔE_def_A.x, ΔE_def_A.y, ΔE_def_A.σ, "Energiile de excitație ale fragmentelor de fisiune", L"\mathrm{\Delta{}E_{def}}", :orangered, :circle, 7, "A", "E [MeV]", 0.0, 2.5, first(ΔE_def_A.x)-0.5, last(ΔE_def_A.x)+0.5);
Plot_E_A = Grafic_scatter(E_sciz_A.x, E_sciz_A.y, E_sciz_A.σ, L"\mathrm{E_{sciziune}}", :navy, :pentagon, 7, Plot_E_A);
Plot_E_A = Grafic_scatter(E_excitatie_A.x, E_excitatie_A.y, E_excitatie_A.σ, L"\mathrm{E}^*", :magenta, :dtriangle, 7, Plot_E_A);
Plot_E_A = Grafic_unire_linie(ΔE_def_A.x, ΔE_def_A.y, ΔE_def_A.σ, Plot_E_A, "", :orangered);
Plot_E_A = Grafic_unire_linie(E_sciz_A.x, E_sciz_A.y, E_sciz_A.σ, Plot_E_A, "", :navy);
Plot_E_A = Grafic_unire_linie(E_excitatie_A.x, E_excitatie_A.y, E_excitatie_A.σ, Plot_E_A, "", :magenta);
Grafic_afisare(Plot_E_A, "Energii_excitatie_A");

Plot_R_A = Grafic_scatter(R_A.x, R_A.y, R_A.σ, latexstring("Raportul de partiționare R = \$\\mathrm{E^*_H/TXE}\$"), "", :orangered, :circle, 7, L"\mathrm{A_H}", "R(A)", 0.9, 1.05, first(R_A.x)-0.5, last(R_A.x)+0.5);
Plot_R_A = Grafic_unire_linie(R_A.x, R_A.y, R_A.σ, Plot_R_A, "", :orangered);
Plot_R_A = Grafic_linie_medie_orizontal(Plot_R_A, 0.5, "", :olive);
Grafic_afisare(Plot_R_A, "R_A");

ν_pereche_med = Medie_distributie_yield(ν_pereche_A, y_A, firstindex(ν_pereche_A.x), lastindex(ν_pereche_A.x));
Plot_ν_pereche_A = Grafic_scatter(ν_pereche_A.x, ν_pereche_A.y, ν_pereche_A.σ, "Multiplicitatea neutronică a perechilor de fragmente în emisia promptă", "Modelare la sciziune", :red, :circle, 10, L"\mathrm{A_H}", latexstring("\$\\nu_{\\mathrm{pereche}}\$"), 0.6, 1.05, first(ν_pereche_A.x)-0.5, last(ν_pereche_A.x)+0.5);
Plot_ν_pereche_A = Grafic_scatter(ν_pereche_Gook.x, ν_pereche_Gook.y, ν_pereche_Gook.σ, "Gook", :blue, :xcross, 6.5, Plot_ν_pereche_A);
Plot_ν_pereche_A = Grafic_scatter(ν_pereche_Maslin.x, ν_pereche_Maslin.y, ν_pereche_Maslin.σ, "Maslin", :green, :star5, 6.5, Plot_ν_pereche_A);
Plot_ν_pereche_A = Grafic_scatter(ν_pereche_Nishio.x, ν_pereche_Nishio.y, ν_pereche_Nishio.σ, "Nishio", :purple, :dtriangle, 6.5, Plot_ν_pereche_A);
Plot_ν_pereche_A = Grafic_scatter(ν_pereche_Vorobyev.x, ν_pereche_Vorobyev.y, ν_pereche_Vorobyev.σ, "Vorobyev", :gold, :diamond, 6.5, Plot_ν_pereche_A);
Plot_ν_pereche_A = Grafic_unire_linie(ν_pereche_A.x, ν_pereche_A.y, ν_pereche_A.σ, Plot_ν_pereche_A, "", :red);
Plot_ν_pereche_A = Grafic_textbox(135, 4.5, Plot_ν_pereche_A, latexstring("\$<\\nu_{\\mathrm{pereche}}>\$ = $(round(ν_pereche_med[1], digits = 3)) \$^\\pm\$ $(round(ν_pereche_med[2], digits = 5))"));
Grafic_afisare(Plot_ν_pereche_A, "Multiplicitati_neutronice_pereche");

ν_med = Medie_distributie_yield(ν_A, y_A, firstindex(ν_A.x), lastindex(ν_A.x));
Plot_ν_A = Grafic_scatter(ν_A.x, ν_A.y, ν_A.σ, "Multiplicitatea neutronică a fragmentelor de fisiune în emisia promptă", "Modelare la sciziune", :red, :circle, 10, "A", latexstring("\$\\nu\$"), 0.0, 1.35, first(ν_A.x)-2.5, last(ν_A.x)+2.5);
Plot_ν_A = Grafic_scatter(ν_Gook.x, ν_Gook.y, ν_Gook.σ, "Gook", :blue, :xcross, 5.5, Plot_ν_A);
Plot_ν_A = Grafic_scatter(ν_Maslin.x, ν_Maslin.y, ν_Maslin.σ, "Maslin", :green, :star5, 5.5, Plot_ν_A);
Plot_ν_A = Grafic_scatter(ν_Nishio.x, ν_Nishio.y, ν_Nishio.σ, "Nishio", :purple, :dtriangle, 5.5, Plot_ν_A);
Plot_ν_A = Grafic_scatter(ν_Vorobyev.x, ν_Vorobyev.y, ν_Vorobyev.σ, "Vorobyev", :gold, :diamond, 5.5, Plot_ν_A);
Plot_ν_A = Grafic_unire_linie(ν_A.x, ν_A.y, ν_A.σ, Plot_ν_A, "", :red);
Plot_ν_A = Grafic_textbox(89, 2.8, Plot_ν_A, latexstring("\$<\\nu>\$ = $(round(ν_med[1], digits = 3)) \$\\pm\$ $(round(ν_med[2], digits = 5))"));
Grafic_afisare(Plot_ν_A, "Multiplicitati_neutronice_fragmente");

Plot_rapoarte_A = Grafic_scatter(R_A.x, R_A.y, R_A.σ, latexstring("Raportul \$\\mathrm{E^*_H/TXE}\$ comparat cu rapoartele \$\\nu_\\mathrm{H}/\\nu_\\mathrm{pereche}\$"), latexstring("R = \$\\mathrm{E^*_H/TXE}\$ - Modelare la sciziune"), :red, :utriangle, 8, L"\mathrm{A_H}", "Rapoarte de partiționare", 0.0, 1.6, first(R_A.x)-0.5, last(R_A.x)+0.5);
Plot_rapoarte_A = Grafic_scatter(raport_νH_νPair_Model.x, raport_νH_νPair_Model.y, raport_νH_νPair_Model.σ, latexstring("\$\\nu_\\mathrm{H}/\\nu_\\mathrm{pereche}\$ - Modelare la sciziune"), :blue, :dtriangle, 8, Plot_rapoarte_A);
Plot_rapoarte_A = Grafic_scatter(raport_νH_νPair_Gook.x, raport_νH_νPair_Gook.y, raport_νH_νPair_Gook.σ, latexstring("\$\\nu_\\mathrm{H}/\\nu_\\mathrm{pereche}\$ - Gook"), :gold, :xcross, 5, Plot_rapoarte_A);
Plot_rapoarte_A = Grafic_scatter(raport_νH_νPair_Nishio.x, raport_νH_νPair_Nishio.y, raport_νH_νPair_Nishio.σ, latexstring("\$\\nu_\\mathrm{H}/\\nu_\\mathrm{pereche}\$ - Nishio"), :green, :star5, 5, Plot_rapoarte_A);
Plot_rapoarte_A = Grafic_scatter(raport_νH_νPair_Maslin.x, raport_νH_νPair_Maslin.y, raport_νH_νPair_Maslin.σ, latexstring("\$\\nu_\\mathrm{H}/\\nu_\\mathrm{pereche}\$ - Maslin"), :purple, :circle, 5, Plot_rapoarte_A);
Plot_rapoarte_A = Grafic_scatter(raport_νH_νPair_Vorobyev.x, raport_νH_νPair_Vorobyev.y, raport_νH_νPair_Vorobyev.σ, latexstring("\$\\nu_\\mathrm{H}/\\nu_\\mathrm{pereche}\$ - Vorobyev"), :orange, :diamond, 5, Plot_rapoarte_A);
Plot_rapoarte_A = Grafic_linie_medie_orizontal(Plot_rapoarte_A, 0.5, "", :olive);
Grafic_afisare(Plot_rapoarte_A, "Rapoarte_A");