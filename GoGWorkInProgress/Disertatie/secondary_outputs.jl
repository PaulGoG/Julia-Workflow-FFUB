#=
This part of the program processes the main output so that
it can be compared with experimental data (averaging main output data over experimental yield distributions)
=#
#####
#Compute Y(A,Z,TKE) form experimental Y(A,TKE) data
function Process_yield_data(A_0, Z_0, fragmdomain::Distribution, dY::DataFrame)
    y_A_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    #Completes fragmentation domain in case data provided only for HF
    if !isassigned(dY.A[dY.A .< A_0/2], 1)
        for A in sort(unique(dY.A), rev = true)
            for Z in fragmdomain.Z[fragmdomain.A .== A_0 - A]
                P_A_Z = fragmdomain.Value[(fragmdomain.A .== A_0 - A) .& (fragmdomain.Z .== Z)][1]
                for TKE in sort(unique(dY.TKE[dY.A .== A]))
                    val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    push!(y_A_Z_TKE.A, A_0 - A)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Value, val)
                    push!(y_A_Z_TKE.σ, σ)
                end
            end
        end
    end
    for A in sort(unique(dY.A))
        for Z in fragmdomain.Z[fragmdomain.A .== A]
            P_A_Z = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
            for TKE in sort(unique(dY.TKE[dY.A .== A]))
                if !isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                    val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    push!(y_A_Z_TKE.A, A)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Value, val)
                    push!(y_A_Z_TKE.σ, σ)
                end
            end
        end
    end
    #Renormalizes data
    f = 200/sum(y_A_Z_TKE.Value)
    y_A_Z_TKE.Value .= y_A_Z_TKE.Value .* f
    y_A_Z_TKE.σ .= y_A_Z_TKE.σ .* f
    return y_A_Z_TKE
end
#Average q(A,Z,TKE) over p(Z,A) so it becomes q(A,TKE)
function Average_over_Z(q_A_Z_TKE, fragmdomain::Distribution)
    q_A_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in sort(unique(q_A_Z_TKE.A))
        for TKE in unique(q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A)])
            Denominator = 0.0
            Numerator = 0.0
            for Z in q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    P_A_Z = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
                    Denominator += P_A_Z
                    Numerator += P_A_Z *value
                end
            end
            if Denominator > 0
                push!(q_A_TKE.A, A)
                push!(q_A_TKE.TKE, TKE)
                push!(q_A_TKE.Value, Numerator/Denominator)
            end
        end
    end
    return q_A_TKE
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) so it becomes q(A)
function Average_over_TKE_Z(q_A_Z_TKE, y_A_Z_TKE::Distribution)
    q_A = Distribution_unidym(Int[], Float64[], Float64[])
    for A in sort(unique(q_A_Z_TKE.A))
        Denominator = 0.0
        Numerator = 0.0
        for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A)])
            for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE *value
                    end
                end
            end
        end
        if Denominator > 0
            push!(q_A.Argument, A)
            push!(q_A.Value, Numerator/Denominator)
        end
    end
    return q_A
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) so it becomes q(TKE)
function Average_over_A_Z(q_A_Z_TKE, y_A_Z_TKE::Distribution)
    q_TKE = Distribution_unidym(Float64[], Float64[], Float64[])
    for TKE in sort(unique(q_A_Z_TKE.TKE))
        Denominator = 0.0
        Numerator = 0.0
        for A in unique(q_A_Z_TKE.A[(q_A_Z_TKE.TKE .== TKE)])
            for Z in q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE *value
                    end
                end
            end
        end
        if Denominator > 0
            push!(q_TKE.Argument, TKE)
            push!(q_TKE.Value, Numerator/Denominator)
        end
    end
    return q_TKE
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) to get average value <q>
function Average_value(q_A_Z_TKE, y_A_Z_TKE::Distribution, mass_number_range)
    Denominator = 0.0
    Numerator = 0.0
    for A in mass_number_range
        for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A)])
            for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE *value
                    end
                end
            end
        end
    end
    return Numerator/Denominator
end
#Average q(Argument) over Y(Argument) to get average value <q>
function Average_value(q, y::Distribution_unidym, Argument_range)
    Denominator = 0.0
    Numerator = 0.0
    for Argument in Argument_range
        if isassigned(q.Value[q.Argument .== Argument], 1) && isassigned(y.Value[y.Argument .== Argument], 1)
            value = q.Value[q.Argument .== Argument][1]
            Y = y.Value[y.Argument .== Argument][1]
            Denominator += Y
            Numerator += Y *value
        end
    end
    return Numerator/Denominator
end
#Get L-H pair value from given q(A) distribution and A_H 
function Pair_value(q_A::Distribution_unidym, A_0, A_H)
    val_L = q_A.Value[q_A.Argument .== A_0 - A_H][1]
    val_H = q_A.Value[q_A.Argument .== A_H][1]
    return val_L + val_H
end
#Obtain Y(Z,Aₚ,TKE), Y(Z,Aₚ) & Y(Aₚ) distributions from Y(A,Z,TKE) & n(A,Z,TKE)
function Yield_post_neutron(y_A_Z_TKE::Distribution, n_A_Z_TKE)
    y_Aₚ_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    y_Aₚ_Z = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in unique(y_A_Z_TKE.A)
        for Z in unique(y_A_Z_TKE.Z[y_A_Z_TKE.A .== A])
            for TKE in y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z)]
                Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                σY_A_Z_TKE = y_A_Z_TKE.σ[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                if isassigned(n_A_Z_TKE.Value[(n_A_Z_TKE.A .== A) .& (n_A_Z_TKE.Z .== Z) .& (n_A_Z_TKE.TKE .== TKE)], 1)
                    n = n_A_Z_TKE.Value[(n_A_Z_TKE.A .== A) .& (n_A_Z_TKE.Z .== Z) .& (n_A_Z_TKE.TKE .== TKE)][1]
                    Aₚ = A - n
                else 
                    Aₚ = A 
                end
                if !isassigned(y_Aₚ_Z_TKE.Value[(y_Aₚ_Z_TKE.A .== Aₚ) .& (y_Aₚ_Z_TKE.Z .== Z) .& (y_Aₚ_Z_TKE.TKE .== TKE)], 1)
                    push!(y_Aₚ_Z_TKE.A, Aₚ)
                    push!(y_Aₚ_Z_TKE.Z, Z)
                    push!(y_Aₚ_Z_TKE.TKE, TKE)
                    push!(y_Aₚ_Z_TKE.Value, Y_A_Z_TKE)
                    push!(y_Aₚ_Z_TKE.σ, σY_A_Z_TKE)
                else
                    y_Aₚ_Z_TKE.Value[(y_Aₚ_Z_TKE.A .== Aₚ) .& (y_Aₚ_Z_TKE.Z .== Z) .& (y_Aₚ_Z_TKE.TKE .== TKE)] .+= Y_A_Z_TKE
                    y_Aₚ_Z_TKE.σ[(y_Aₚ_Z_TKE.A .== Aₚ) .& (y_Aₚ_Z_TKE.Z .== Z) .& (y_Aₚ_Z_TKE.TKE .== TKE)] .= sqrt(sum(y_Aₚ_Z_TKE.σ[(y_Aₚ_Z_TKE.A .== Aₚ) .& (y_Aₚ_Z_TKE.Z .== Z) .& (y_Aₚ_Z_TKE.TKE .== TKE)].^2) + σY_A_Z_TKE^2)
                end
            end
        end
    end
    for A in sort(unique(y_Aₚ_Z_TKE.A))
        for Z in sort(unique(y_Aₚ_Z_TKE.Z[(y_Aₚ_Z_TKE.A .== A)]))
            Y_Aₚ_Z = sum(y_Aₚ_Z_TKE.Value[(y_Aₚ_Z_TKE.A .== A) .& (y_Aₚ_Z_TKE.Z .== Z)])
            σY_Aₚ_Z = sqrt(sum(y_Aₚ_Z_TKE.σ[(y_Aₚ_Z_TKE.A .== A) .& (y_Aₚ_Z_TKE.Z .== Z)].^2))
            push!(y_Aₚ_Z.A, A)
            push!(y_Aₚ_Z.Z, Z)
            push!(y_Aₚ_Z.Value, Y_Aₚ_Z)
            push!(y_Aₚ_Z.σ, σY_Aₚ_Z)
        end
    end
    return y_Aₚ_Z_TKE, y_Aₚ_Z
end
#Obtain vectorized singular distributions Y(A), Y(N), Y(Z), Y(TKE), TKE(AH) & KE(A) from Y(A,Z,TKE)
function Singular_yield_distributions(y_A_Z_TKE::Distribution, A_0, A_H_min)
    y_A = Distribution_unidym(Int[], Float64[], Float64[])
    y_Z = Distribution_unidym(Int[], Float64[], Float64[])
    y_N = Distribution_unidym(Int[], Float64[], Float64[])
    y_TKE = Distribution_unidym(Float64[], Float64[], Float64[])
    tke_AH = Distribution_unidym(Int[], Float64[], Float64[])
    ke_A = Distribution_unidym(Int[], Float64[], Float64[])
    #Y(A)
    for A in sort(unique(y_A_Z_TKE.A))
        Y_A = sum(y_A_Z_TKE.Value[y_A_Z_TKE.A .== A])
        σY_A = sqrt(sum(y_A_Z_TKE.σ[y_A_Z_TKE.A .== A] .^2))
        push!(y_A.Argument, A)
        push!(y_A.Value, Y_A)
        push!(y_A.σ, σY_A)
    end
    #Y(Z)
    for Z in sort(unique(y_A_Z_TKE.Z))
        Y_Z = sum(y_A_Z_TKE.Value[y_A_Z_TKE.Z .== Z])
        σY_Z = sqrt(sum(y_A_Z_TKE.σ[y_A_Z_TKE.Z .== Z] .^2))
        push!(y_Z.Argument, Z)
        push!(y_Z.Value, Y_Z)
        push!(y_Z.σ, σY_Z)
    end
    #Y(N)
    for A in sort(unique(y_A_Z_TKE.A))
        for Z in unique(y_A_Z_TKE.Z[y_A_Z_TKE.A .== A])
            N = A - Z
            if !isassigned(y_N.Value[y_N.Argument .== N], 1)
                Y_N = sum(y_A_Z_TKE.Value[y_A_Z_TKE.A .- y_A_Z_TKE.Z .== N])
                σY_N = sqrt(sum(y_A_Z_TKE.σ[y_A_Z_TKE.A .- y_A_Z_TKE.Z .== N].^2))
                push!(y_N.Argument, N)
                push!(y_N.Value, Y_N)
                push!(y_N.σ, σY_N)
            end
        end
    end
    #Y(TKE)
    for TKE in sort(unique(y_A_Z_TKE.TKE))
        Y_TKE = sum(y_A_Z_TKE.Value[y_A_Z_TKE.TKE .== TKE])
        σY_TKE = sqrt(sum(y_A_Z_TKE.σ[y_A_Z_TKE.TKE .== TKE] .^2))
        push!(y_TKE.Argument, TKE)
        push!(y_TKE.Value, Y_TKE)
        push!(y_TKE.σ, σY_TKE)
    end
    #TKE(AH)
    for A_H in sort(unique(y_A_Z_TKE.A[y_A_Z_TKE.A .>= A_H_min]))
        Numerator = 0.0
        𝚺_σ² = 0.0
        Denominator = sum(y_A_Z_TKE.Value[y_A_Z_TKE.A .== A_H])
        for TKE in sort(unique(y_A_Z_TKE.TKE[y_A_Z_TKE.A .== A_H]))
            Y_A_TKE = sum(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.TKE .== TKE)])
            Numerator += TKE *Y_A_TKE
        end
        TKE_AH = Numerator/Denominator
        for TKE in sort(unique(y_A_Z_TKE.TKE[y_A_Z_TKE.A .== A_H]))
            σY_A_TKE = sqrt(sum(y_A_Z_TKE.σ[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.TKE .== TKE)] .^2))
            𝚺_σ² += (TKE - TKE_AH)^2 *σY_A_TKE^2
        end
        push!(tke_AH.Argument, A_H)
        push!(tke_AH.Value, TKE_AH)
        push!(tke_AH.σ, sqrt(𝚺_σ²)/Denominator)
    end  
    #KE(A)
    for A_H in sort(unique(tke_AH.Argument), rev = true)
        TKE_A = tke_AH.Value[tke_AH.Argument .== A_H][1]
        σTKE_A = tke_AH.σ[tke_AH.Argument .== A_H][1]
        KE_A = TKE_A *A_H /A_0
        push!(ke_A.Argument, A_0 - A_H)
        push!(ke_A.Value, KE_A)
        push!(ke_A.σ, KE_A *σTKE_A /TKE_A)
    end
    for A_H in unique(tke_AH.Argument)
        if A_H != A_0 - A_H
            TKE_A = tke_AH.Value[tke_AH.Argument .== A_H][1]
            σTKE_A = tke_AH.σ[tke_AH.Argument .== A_H][1]
            KE_A = TKE_A *(A_0 - A_H) /A_0
            push!(ke_A.Argument, A_H)
            push!(ke_A.Value, KE_A)
            push!(ke_A.σ, KE_A *σTKE_A /TKE_A)
        end
    end
    Sort_q_Argument(y_A)
    Sort_q_Argument(y_Z)
    Sort_q_Argument(y_N)
    Sort_q_Argument(y_TKE)
    Sort_q_Argument(tke_AH)
    Sort_q_Argument(ke_A)
    return y_A, y_Z, y_N, y_TKE, tke_AH, ke_A
end
#Get <Argument> from Yield(Argument) singular yield distribution
function Average_yield_argument(yield::Distribution_unidym, argument_range)
    Numerator = 0.0
    Denominator = 0.0
    𝚺_σ² = 0.0
    for Argument in argument_range
        Y = yield.Value[yield.Argument .== Argument][1]
        Numerator += Argument *Y
        Denominator += Y
    end
    Average_arg = Numerator/Denominator
    for Argument in argument_range
        σY = yield.σ[yield.Argument .== Argument][1]
        𝚺_σ² += (Argument - Average_arg)^2 *σY^2
    end    
    return Average_arg, sqrt(𝚺_σ²)/Denominator
end
#Obtain vectorized singular distributions Q(AH), TXE(AH)
function Vectorized_TXE_Q_AH(A_0, Z_0, fission_type::String, E_incident, y_A_Z_TKE::Distribution, A_H_range, dm::DataFrame)
    Q_AH = Distribution_unidym(Int[], Float64[], Float64[])
    txe_AH = Distribution_unidym(Int[], Float64[], Float64[])
    E_CN = Compound_nucleus_energy(fission_type, A_0, Z_0, E_incident, dm)
    for A_H in A_H_range
        Denominator_Q = 0.0
        Numerator_Q = 0.0
        Denominator_TXE = 0.0
        Numerator_TXE = 0.0
        for Z_H in unique(y_A_Z_TKE.Z[(y_A_Z_TKE.A .== A_H)])
            Q_A_Z = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
            for TKE in y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H)]
                TXE_A_Z_TKE = Total_excitation_energy(Q_A_Z[1], Q_A_Z[2], TKE, 0.0, E_CN[1], E_CN[2])
                if !isnan(TXE_A_Z_TKE[1])
                    Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H) .& (y_A_Z_TKE.TKE .== TKE)][1]
                    Numerator_TXE += Y_A_Z_TKE *TXE_A_Z_TKE[1]
                    Denominator_TXE += Y_A_Z_TKE
                end
            end
            Y_A_Z = sum(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H)])
            Numerator_Q += Y_A_Z *Q_A_Z[1]
            Denominator_Q += Y_A_Z
        end
        if Denominator_TXE > 0
            push!(txe_AH.Argument, A_H)
            push!(txe_AH.Value, Numerator_TXE/Denominator_TXE)
        end
        if Denominator_Q > 0
            push!(Q_AH.Argument, A_H)
            push!(Q_AH.Value, Numerator_Q/Denominator_Q)
        end
    end
    return Q_AH, txe_AH
end
#Compute probabilities of occurance P(q) of any vectorized generic quantity q(A,Z,TKE) with value step of Δq
function Probability_of_occurrence(q_A_Z_TKE, y_A_Z_TKE::Distribution, Δq::Number)
    P = Distribution_unidym(Float64[], Float64[], Float64[])
    q_values = filter(!isnan, q_A_Z_TKE.Value)
    if isassigned(q_values, 1)
        lower_q = floor(minimum(q_values))
        upper_q = lower_q + Δq
        while lower_q <= maximum(q_values)
            Frequency = 0.0
            for A in unique(q_A_Z_TKE.A[(q_A_Z_TKE.Value .>= lower_q) .& (q_A_Z_TKE.Value .< upper_q)])
                for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Value .>= lower_q) .& (q_A_Z_TKE.Value .< upper_q)])
                    for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z) .& (q_A_Z_TKE.Value .>= lower_q) .& (q_A_Z_TKE.Value .< upper_q)]
                        if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                            Frequency += y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        end
                    end
                end
            end
            push!(P.Argument, lower_q)
            push!(P.Value, Frequency)
            lower_q = upper_q
            upper_q += Δq
        end
        f = 100/sum(P.Value)
        P.Value .*= f
    else
        push!(P.Argument, NaN)
        push!(P.Value, NaN)
    end
    return P
end
#Sort q(Argument) by ascending argument for plotting
function Sort_q_Argument(q)
    aux_Argument = sort(unique(q.Argument))
    aux_Value = [first(q.Value[q.Argument .== Argument]) for Argument in aux_Argument]
    for index in eachindex(aux_Argument)
        q.Argument[index] = aux_Argument[index]
        q.Value[index] = aux_Value[index]
    end
    if isassigned(filter(!isnan, q.σ), 1)
        aux_σ = [first(q.σ[q.Argument .== Argument]) for Argument in aux_Argument]
        for index in eachindex(aux_Argument)
            q.σ[index] = aux_σ[index]
        end
    end
end
#####
println("*averaging data over $yield_distribution_filename experimental Yield distribution") 

y_A_Z_TKE = Process_yield_data(A₀, Z₀, fragmdomain, Yield_data)

if secondary_output_Yield == "YES"
    y_A, y_Z, y_N, y_TKE, tke_AH, ke_A = Singular_yield_distributions(y_A_Z_TKE, A₀, A_H_min)
    if !isdir("output_data/Yield/")
        mkdir("output_data/Yield/")
    end
    if isassigned(filter(!isnan, y_A_Z_TKE.σ), 1)
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_A_Z_TKE.OUT", 
            DataFrame(A = y_A_Z_TKE.A, Z = y_A_Z_TKE.Z, TKE = y_A_Z_TKE.TKE, Y = y_A_Z_TKE.Value, σ = y_A_Z_TKE.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_A.OUT", 
            DataFrame(A = y_A.Argument, Y = y_A.Value, σ = y_A.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_Z.OUT", 
            DataFrame(Z = y_Z.Argument, Y = y_Z.Value, σ = y_Z.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_N.OUT", 
            DataFrame(N = y_N.Argument, Y = y_N.Value, σ = y_N.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_TKE.OUT", 
            DataFrame(TKE = y_TKE.Argument, Y = y_TKE.Value, σ = y_TKE.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_TKE_AH.OUT", 
            DataFrame(A_H = tke_AH.Argument, TKE = tke_AH.Value, σ = tke_AH.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_KE_A.OUT", 
            DataFrame(A = ke_A.Argument, KE = ke_A.Value, σ = ke_A.σ), 
            writeheader=true, newline="\r\n", delim=' '
        )
    else
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_A_Z_TKE.OUT", 
            DataFrame(A = y_A_Z_TKE.A, Z = y_A_Z_TKE.Z, TKE = y_A_Z_TKE.TKE, Y = y_A_Z_TKE.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_A.OUT", 
            DataFrame(A = y_A.Argument, Y = y_A.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_Z.OUT", 
            DataFrame(Z = y_Z.Argument, Y = y_Z.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_N.OUT", 
            DataFrame(N = y_N.Argument, Y = y_N.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_Y_TKE.OUT", 
            DataFrame(TKE = y_TKE.Argument, Y = y_TKE.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_TKE_AH.OUT", 
            DataFrame(A_H = tke_AH.Argument, TKE = tke_AH.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
        CSV.write(
            "output_data/Yield/$(fissionant_nucleus_identifier)_KE_A.OUT", 
            DataFrame(A = ke_A.Argument, KE = ke_A.Value), 
            writeheader=true, newline="\r\n", delim=' '
        )
    end
end
if secondary_output_ν == "YES"
    ν_A_Z_TKE = Neutron_multiplicity_A_Z_TKE(DataFrame(
        A = Raw_output_datafile.A,
        Z = Raw_output_datafile.Z,
        TKE = Raw_output_datafile.TKE,
        No_Sequence = Raw_output_datafile.No_Sequence
        )
    )
    if !isdir("output_data/nu/")
        mkdir("output_data/nu/")
    end
    max_seq_A_Z_TKE = Maximum_sequences_A_Z_TKE(DataFrame(
        A = Raw_output_datafile.A,
        Z = Raw_output_datafile.Z,
        TKE = Raw_output_datafile.TKE,
        No_Sequence = Raw_output_datafile.No_Sequence
        )
    )
    ν_A_TKE = Average_over_Z(ν_A_Z_TKE, fragmdomain)
    CSV.write(
        "output_data/nu/$(fissionant_nucleus_identifier)_nu_A_TKE.OUT", 
        DataFrame(A = ν_A_TKE.A, TKE = ν_A_TKE.TKE, ν = ν_A_TKE.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    ν_A = Average_over_TKE_Z(ν_A_Z_TKE, y_A_Z_TKE)
    CSV.write(
        "output_data/nu/$(fissionant_nucleus_identifier)_nu_A.OUT", 
        DataFrame(A = ν_A.Argument, ν = ν_A.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    ν_AH_Pair = DataFrame(
        Argument =  ν_A.Argument[ν_A.Argument .>= A_H_min],
        Value = [Pair_value(ν_A, A₀, A_H) for A_H in ν_A.Argument[ν_A.Argument .>= A_H_min]]
    )
    CSV.write(
        "output_data/nu/$(fissionant_nucleus_identifier)_nu_AH_Pair.OUT", 
        DataFrame(A = ν_AH_Pair.Argument, ν_Pair = ν_AH_Pair.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    ν_TKE = Average_over_A_Z(ν_A_Z_TKE, y_A_Z_TKE)
    CSV.write(
        "output_data/nu/$(fissionant_nucleus_identifier)_nu_Pair_TKE.OUT", 
        DataFrame(TKE = ν_TKE.Argument, ν = ν_TKE.Value .* 2), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_ν = Probability_of_occurrence(ν_A_Z_TKE, y_A_Z_TKE, 1)
    CSV.write(
        "output_data/nu/$(fissionant_nucleus_identifier)_P_nu.OUT", 
        DataFrame(ν = probability_ν.Argument, P = probability_ν.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    #=
    probability_ν_L = Probability_of_occurrence(
        DataFrame(
            A = ν_A_Z_TKE.A[ν_A_Z_TKE.A .<= A_H_min], 
            Z = ν_A_Z_TKE.Z[ν_A_Z_TKE.A .<= A_H_min], 
            TKE = ν_A_Z_TKE.TKE[ν_A_Z_TKE.A .<= A_H_min],
            Value = ν_A_Z_TKE.Value[ν_A_Z_TKE.A .<= A_H_min]
        ), y_A_Z_TKE, 1
    )
    CSV.write(
        "output_data/nu/$(fissionant_nucleus_identifier)_P_nu_LF.OUT", 
        DataFrame(ν = probability_ν_L.Argument, P = probability_ν_L.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_ν_H = Probability_of_occurrence(
        DataFrame(
            A = ν_A_Z_TKE.A[ν_A_Z_TKE.A .>= A_H_min], 
            Z = ν_A_Z_TKE.Z[ν_A_Z_TKE.A .>= A_H_min], 
            TKE = ν_A_Z_TKE.TKE[ν_A_Z_TKE.A .>= A_H_min],
            Value = ν_A_Z_TKE.Value[ν_A_Z_TKE.A .>= A_H_min]
        ), y_A_Z_TKE, 1
    )
    CSV.write(
        "output_data/nu/$(fissionant_nucleus_identifier)_P_nu_HF.OUT", 
        DataFrame(ν = probability_ν_H.Argument, P = probability_ν_H.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    =#
    if secondary_output_Ap == "YES"
        y_Ap_Z_TKE, y_Ap_Z = Yield_post_neutron(y_A_Z_TKE, max_seq_A_Z_TKE)
        Ap_H_min = A_H_min - ceil(ν_A.Value[ν_A.Argument .== A_H_min][1])
        y_Ap, y_Zp, y_Np, y_TKEp, tke_AHp, ke_Ap = Singular_yield_distributions(y_Ap_Z_TKE, A₀, Ap_H_min)
        if !isdir("output_data/Yield_Ap/")
            mkdir("output_data/Yield_Ap/")
        end
        if !isdir("output_data/Yield_Ap/Yield_Ap_Z/")
            mkdir("output_data/Yield_Ap/Yield_Ap_Z/")
        end
        if isassigned(filter(!isnan, y_Ap_Z_TKE.σ), 1)
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_Y_Ap.OUT", 
                DataFrame(Aₚ = y_Ap.Argument, Y = y_Ap.Value, σ = y_Ap.σ), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_Y_Ap_Z.OUT", 
                DataFrame(Aₚ = y_Ap_Z.A, Z = y_Ap_Z.Z, Y = y_Ap_Z.Value, σ = y_Ap_Z.σ), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_Y_Np.OUT", 
                DataFrame(N = y_Np.Argument, Y = y_Np.Value, σ = y_Np.σ), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_TKE_AHp.OUT", 
                DataFrame(A_H = tke_AHp.Argument, TKE = tke_AHp.Value, σ = tke_AHp.σ), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_KE_Ap.OUT", 
                DataFrame(A = ke_Ap.Argument, KE = ke_Ap.Value, σ = ke_Ap.σ), 
                writeheader=true, newline="\r\n", delim=' '
            )
            for Z in unique(y_Ap_Z.Z)
                CSV.write(
                    "output_data/Yield_Ap/Yield_Ap_Z/$(fissionant_nucleus_identifier)_Y_Ap_$(Z).OUT", 
                    DataFrame(Aₚ = y_Ap_Z.A[y_Ap_Z.Z .== Z], Y = y_Ap_Z.Value[y_Ap_Z.Z .== Z], σ = y_Ap_Z.σ[y_Ap_Z.Z .== Z]), 
                    writeheader=true, newline="\r\n", delim=' '
                )
            end
        else
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_Y_Ap.OUT", 
                DataFrame(Aₚ = y_Ap.Argument, Y = y_Ap.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_Y_Ap_Z.OUT", 
                DataFrame(Aₚ = y_Ap_Z.A, Z = y_Ap_Z.Z, Y = y_Ap_Z.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_Y_Np.OUT", 
                DataFrame(N = y_Np.Argument, Y = y_Np.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_TKE_AHp.OUT", 
                DataFrame(A_H = tke_AHp.Argument, TKE = tke_AHp.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            CSV.write(
                "output_data/Yield_Ap/$(fissionant_nucleus_identifier)_KE_Ap.OUT", 
                DataFrame(A = ke_Ap.Argument, KE = ke_Ap.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            for Z in sort(unique(y_Ap_Z.Z))
                CSV.write(
                    "output_data/Yield_Ap/Yield_Ap_Z/$(fissionant_nucleus_identifier)_Y_Ap_$(Z).OUT", 
                    DataFrame(Aₚ = y_Ap_Z.A[y_Ap_Z.Z .== Z], Y = y_Ap_Z.Value[y_Ap_Z.Z .== Z]), 
                    writeheader=true, newline="\r\n", delim=' '
                )
            end
        end
    end
    if secondary_output_Tₖ == "YES"
        if !isdir("output_data/P_T_k/")
            mkdir("output_data/P_T_k/")
        end
        if !isdir("output_data/T_k_A/")
            mkdir("output_data/T_k_A/")
        end
        if secondary_output_Eᵣ == "YES"
            if !isdir("output_data/P_Er/")
                mkdir("output_data/P_Er/")
            end
        end
        for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
            Tₖ_A_Z_TKE = DataFrame(
                A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                Value = Raw_output_datafile.Tₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
            )
            Tₖ_A = Average_over_TKE_Z(Tₖ_A_Z_TKE, y_A_Z_TKE)
            CSV.write(
                "output_data/T_k_A/$(fissionant_nucleus_identifier)_T_$(k)_A.OUT", 
                DataFrame(A = Tₖ_A.Argument, T = Tₖ_A.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_Tₖ = Probability_of_occurrence(Tₖ_A_Z_TKE, y_A_Z_TKE, ΔTₖ)
            CSV.write(
                "output_data/P_T_k/$(fissionant_nucleus_identifier)_P_T_$(k).OUT", 
                DataFrame(Tₖ = probability_Tₖ.Argument, P = probability_Tₖ.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_Tₖ_L = Probability_of_occurrence(
                DataFrame(
                    A = Tₖ_A_Z_TKE.A[Tₖ_A_Z_TKE.A .<= A_H_min], 
                    Z = Tₖ_A_Z_TKE.Z[Tₖ_A_Z_TKE.A .<= A_H_min], 
                    TKE = Tₖ_A_Z_TKE.TKE[Tₖ_A_Z_TKE.A .<= A_H_min],
                    Value = Tₖ_A_Z_TKE.Value[Tₖ_A_Z_TKE.A .<= A_H_min]
                ), y_A_Z_TKE, ΔTₖ
                )
            CSV.write(
                "output_data/P_T_k/$(fissionant_nucleus_identifier)_P_T_$(k)_LF.OUT", 
                DataFrame(Tₖ = probability_Tₖ_L.Argument, P = probability_Tₖ_L.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_Tₖ_H = Probability_of_occurrence(
                DataFrame(
                    A = Tₖ_A_Z_TKE.A[Tₖ_A_Z_TKE.A .>= A_H_min], 
                    Z = Tₖ_A_Z_TKE.Z[Tₖ_A_Z_TKE.A .>= A_H_min], 
                    TKE = Tₖ_A_Z_TKE.TKE[Tₖ_A_Z_TKE.A .>= A_H_min],
                    Value = Tₖ_A_Z_TKE.Value[Tₖ_A_Z_TKE.A .>= A_H_min]
                ), y_A_Z_TKE, ΔTₖ
                )
            CSV.write(
                "output_data/P_T_k/$(fissionant_nucleus_identifier)_P_T_$(k)_HF.OUT", 
                DataFrame(Tₖ = probability_Tₖ_H.Argument, P = probability_Tₖ_H.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            if secondary_output_Eᵣ == "YES"
                aₖ = copy(Raw_output_datafile.aₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)])
                Eᵣ_A_Z_TKE = copy(Tₖ_A_Z_TKE)
                Eᵣ_A_Z_TKE.Value .= Energy_FermiGas.(aₖ, Eᵣ_A_Z_TKE.Value)
                probability_Eᵣ = Probability_of_occurrence(Eᵣ_A_Z_TKE, y_A_Z_TKE, ΔEᵣ)
                CSV.write(
                    "output_data/P_Er/$(fissionant_nucleus_identifier)_P_Er_$(k).OUT", 
                    DataFrame(Eᵣ = probability_Eᵣ.Argument, P = probability_Eᵣ.Value), 
                    writeheader=true, newline="\r\n", delim=' '
                )
                probability_Eᵣ_L = Probability_of_occurrence(
                    DataFrame(
                        A = Eᵣ_A_Z_TKE.A[Eᵣ_A_Z_TKE.A .<= A_H_min], 
                        Z = Eᵣ_A_Z_TKE.Z[Eᵣ_A_Z_TKE.A .<= A_H_min], 
                        TKE = Eᵣ_A_Z_TKE.TKE[Eᵣ_A_Z_TKE.A .<= A_H_min],
                        Value = Eᵣ_A_Z_TKE.Value[Eᵣ_A_Z_TKE.A .<= A_H_min]
                    ), y_A_Z_TKE, ΔEᵣ
                    )
                CSV.write(
                    "output_data/P_Er/$(fissionant_nucleus_identifier)_P_Er_$(k)_LF.OUT", 
                    DataFrame(Eᵣ = probability_Eᵣ_L.Argument, P = probability_Eᵣ_L.Value), 
                    writeheader=true, newline="\r\n", delim=' '
                )
                probability_Eᵣ_H = Probability_of_occurrence(
                    DataFrame(
                        A = Eᵣ_A_Z_TKE.A[Eᵣ_A_Z_TKE.A .>= A_H_min], 
                        Z = Eᵣ_A_Z_TKE.Z[Eᵣ_A_Z_TKE.A .>= A_H_min], 
                        TKE = Eᵣ_A_Z_TKE.TKE[Eᵣ_A_Z_TKE.A .>= A_H_min],
                        Value = Eᵣ_A_Z_TKE.Value[Eᵣ_A_Z_TKE.A .>= A_H_min]
                    ), y_A_Z_TKE, ΔEᵣ
                    )
                CSV.write(
                    "output_data/P_Er/$(fissionant_nucleus_identifier)_P_Er_$(k)_HF.OUT", 
                    DataFrame(Eᵣ = probability_Eᵣ_H.Argument, P = probability_Eᵣ_H.Value), 
                    writeheader=true, newline="\r\n", delim=' '
                )
            end      
        end
    end
    if secondary_output_avg_εₖ == "YES"
        if !isdir("output_data/P_avgE_k/")
            mkdir("output_data/P_avgE_k/")
        end
        for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
            avg_εₖ_A_Z_TKE = DataFrame(
                A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                Value = Raw_output_datafile.Avg_εₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
            )
            probability_avg_εₖ = Probability_of_occurrence(avg_εₖ_A_Z_TKE, y_A_Z_TKE, Δavg_εₖ)
            CSV.write(
                "output_data/P_avgE_k/$(fissionant_nucleus_identifier)_P_avgE_$(k)_SCM.OUT", 
                DataFrame(Avg_εₖ = probability_avg_εₖ.Argument, P = probability_avg_εₖ.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_avg_εₖ_L = Probability_of_occurrence(
                DataFrame(
                    A = avg_εₖ_A_Z_TKE.A[avg_εₖ_A_Z_TKE.A .<= A_H_min], 
                    Z = avg_εₖ_A_Z_TKE.Z[avg_εₖ_A_Z_TKE.A .<= A_H_min], 
                    TKE = avg_εₖ_A_Z_TKE.TKE[avg_εₖ_A_Z_TKE.A .<= A_H_min],
                    Value = avg_εₖ_A_Z_TKE.Value[avg_εₖ_A_Z_TKE.A .<= A_H_min]
                ), y_A_Z_TKE, Δavg_εₖ
            )
            CSV.write(
                "output_data/P_avgE_k/$(fissionant_nucleus_identifier)_P_avgE_$(k)_SCM_LF.OUT", 
                DataFrame(avg_εₖ = probability_avg_εₖ_L.Argument, P = probability_avg_εₖ_L.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
            probability_avg_εₖ_H = Probability_of_occurrence(
                DataFrame(
                    A = avg_εₖ_A_Z_TKE.A[avg_εₖ_A_Z_TKE.A .>= A_H_min], 
                    Z = avg_εₖ_A_Z_TKE.Z[avg_εₖ_A_Z_TKE.A .>= A_H_min], 
                    TKE = avg_εₖ_A_Z_TKE.TKE[avg_εₖ_A_Z_TKE.A .>= A_H_min],
                    Value = avg_εₖ_A_Z_TKE.Value[avg_εₖ_A_Z_TKE.A .>= A_H_min]
                ), y_A_Z_TKE, Δavg_εₖ
            )
            CSV.write(
                "output_data/P_avgE_k/$(fissionant_nucleus_identifier)_P_avgE_$(k)_SCM_HF.OUT", 
                DataFrame(avg_εₖ = probability_avg_εₖ_H.Argument, P = probability_avg_εₖ_H.Value), 
                writeheader=true, newline="\r\n", delim=' '
            )
        end
    end
end
if secondary_output_T == "YES"
    T_A_Z_TKE = SeqAvg_A_Z_TKE(DataFrame(
        A = Raw_output_datafile.A,
        Z = Raw_output_datafile.Z,
        TKE = Raw_output_datafile.TKE,
        No_Sequence = Raw_output_datafile.No_Sequence,
        Value = Raw_output_datafile.Tₖ
    ))
    T_A = Average_over_TKE_Z(T_A_Z_TKE, y_A_Z_TKE)
    CSV.write(
        "output_data/$(fissionant_nucleus_identifier)_T_A.OUT", 
        DataFrame(A = T_A.Argument, T = T_A.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_T = Probability_of_occurrence(T_A_Z_TKE, y_A_Z_TKE, ΔT)
    CSV.write(
        "output_data/$(fissionant_nucleus_identifier)_P_T.OUT", 
        DataFrame(T = probability_T.Argument, P = probability_T.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_T_L = Probability_of_occurrence(
        DataFrame(
            A = T_A_Z_TKE.A[T_A_Z_TKE.A .<= A_H_min], 
            Z = T_A_Z_TKE.Z[T_A_Z_TKE.A .<= A_H_min], 
            TKE = T_A_Z_TKE.TKE[T_A_Z_TKE.A .<= A_H_min],
            Value = T_A_Z_TKE.Value[T_A_Z_TKE.A .<= A_H_min]
        ), y_A_Z_TKE, ΔT
        )
    CSV.write(
        "output_data/$(fissionant_nucleus_identifier)_P_T_LF.OUT", 
        DataFrame(T = probability_T_L.Argument, P = probability_T_L.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_T_H = Probability_of_occurrence(
        DataFrame(
            A = T_A_Z_TKE.A[T_A_Z_TKE.A .>= A_H_min], 
            Z = T_A_Z_TKE.Z[T_A_Z_TKE.A .>= A_H_min], 
            TKE = T_A_Z_TKE.TKE[T_A_Z_TKE.A .>= A_H_min],
            Value = T_A_Z_TKE.Value[T_A_Z_TKE.A .>= A_H_min]
        ), y_A_Z_TKE, ΔT
        )
    CSV.write(
        "output_data/$(fissionant_nucleus_identifier)_P_T_HF.OUT", 
        DataFrame(T = probability_T_H.Argument, P = probability_T_H.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
end
if secondary_output_avg_ε == "YES"
    avg_ε_A_Z_TKE = SeqAvg_A_Z_TKE(DataFrame(
        A = Raw_output_datafile.A,
        Z = Raw_output_datafile.Z,
        TKE = Raw_output_datafile.TKE,
        No_Sequence = Raw_output_datafile.No_Sequence,
        Value = Raw_output_datafile.Avg_εₖ
    ))
    probability_avg_ε = Probability_of_occurrence(avg_ε_A_Z_TKE, y_A_Z_TKE, Δavg_ε)
    CSV.write(
        "output_data/$(fissionant_nucleus_identifier)_P_avgE_SCM.OUT", 
        DataFrame(avg_ε = probability_avg_ε.Argument, P = probability_avg_ε.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_avg_ε_L = Probability_of_occurrence(
        DataFrame(
            A = avg_ε_A_Z_TKE.A[avg_ε_A_Z_TKE.A .<= A_H_min], 
            Z = avg_ε_A_Z_TKE.Z[avg_ε_A_Z_TKE.A .<= A_H_min], 
            TKE = avg_ε_A_Z_TKE.TKE[avg_ε_A_Z_TKE.A .<= A_H_min],
            Value = avg_ε_A_Z_TKE.Value[avg_ε_A_Z_TKE.A .<= A_H_min]
        ), y_A_Z_TKE, Δavg_ε
        )
    CSV.write(
        "output_data/$(fissionant_nucleus_identifier)_P_avgE_SCM_LF.OUT", 
        DataFrame(avg_ε = probability_avg_ε_L.Argument, P = probability_avg_ε_L.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_avg_ε_H = Probability_of_occurrence(
        DataFrame(
            A = avg_ε_A_Z_TKE.A[avg_ε_A_Z_TKE.A .>= A_H_min], 
            Z = avg_ε_A_Z_TKE.Z[avg_ε_A_Z_TKE.A .>= A_H_min], 
            TKE = avg_ε_A_Z_TKE.TKE[avg_ε_A_Z_TKE.A .>= A_H_min],
            Value = avg_ε_A_Z_TKE.Value[avg_ε_A_Z_TKE.A .>= A_H_min]
        ), y_A_Z_TKE, Δavg_ε
        )
    CSV.write(
        "output_data/$(fissionant_nucleus_identifier)_P_avgE_SCM_HF.OUT", 
        DataFrame(avg_ε = probability_avg_ε_H.Argument, P = probability_avg_ε_H.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
end
if secondary_output_TXE_Q == "YES"
    Q_AH, txe_AH = Vectorized_TXE_Q_AH(A₀, Z₀, fission_type, E_incident, y_A_Z_TKE, A_H_range, dmass_excess)
end
if secondary_output_E_excitation == "YES"
    E_excitation_A = Average_over_TKE_Z(E_excitation, y_A_Z_TKE)
    CSV.write(
        "output_data/$(fissionant_nucleus_identifier)_E_excit_A.OUT", 
        DataFrame(A = E_excitation_A.Argument, E = E_excitation_A.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )
    probability_E_excitation = Probability_of_occurrence(E_excitation, y_A_Z_TKE, 1.0)
    CSV.write(
        "output_data/$(fissionant_nucleus_identifier)_P_E_excit.OUT", 
        DataFrame(E = probability_E_excitation.Argument, P = probability_E_excitation.Value), 
        writeheader=true, newline="\r\n", delim=' '
    )       
end
#Write average quantities to file
open("output_data/$(fissionant_nucleus_identifier)_Average_quantities.OUT", "w") do file
    if secondary_output_Yield == "YES"
        avg_A_L = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .<= A_H_min])
        avg_A_H = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .>= A_H_min])
        if !isnan(avg_A_H[2])
            write(file, "<A>_L = $(avg_A_L[1]) ± $(avg_A_L[2])\n")
            write(file, "<A>_H = $(avg_A_H[1]) ± $(avg_A_H[2])\n")
        else
            write(file, "<A>_L = $(avg_A_L[1])\n")
            write(file, "<A>_H = $(avg_A_H[1])\n")
        end
        avg_TKE = Average_yield_argument(y_TKE, y_TKE.Argument)
        if !isnan(avg_TKE[2])
            write(file, "<TKE> = $(avg_TKE[1]) ± $(avg_TKE[2])\n")
        else
            write(file, "<TKE> = $(avg_TKE[1])\n")
        end
        δₑₒ = (sum(y_Z.Value[iseven.(y_Z.Argument)]) - sum(y_Z.Value[isodd.(y_Z.Argument)]))/sum(y_Z.Value)
        σδₑₒ = (1/sum(y_Z.Value)) * sqrt((1 + δₑₒ)^2 * sum(y_Z.σ .^2) + 2*δₑₒ*(sum(y_Z.σ[isodd.(y_Z.Argument)] .^2) - sum(y_Z.σ[iseven.(y_Z.Argument)].^2)))
        if !isnan(σδₑₒ)
            write(file, "δₑₒ = $(δₑₒ *100) ± $(σδₑₒ *100) %\n\n")
        else
            write(file, "δₑₒ = $(δₑₒ *100) %\n\n")
        end
    end
    if secondary_output_ν == "YES"
        avg_ν_L = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_L_range)
        avg_ν_H = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_H_range)
        avg_ν = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_range)
        avg_ν_Pair = Average_value(ν_AH_Pair, y_A, A_H_range)
        write(file, "<ν>_L = $avg_ν_L\n<ν>_H = $avg_ν_H\n<ν> = $avg_ν\n<ν>_pair = $avg_ν_Pair\n\n")
        if secondary_output_Ap == "YES"
            avg_Ap_L = Average_yield_argument(y_Ap, y_Ap.Argument[y_Ap.Argument .<= Ap_H_min])
            avg_Ap_H = Average_yield_argument(y_Ap, y_Ap.Argument[y_Ap.Argument .>= Ap_H_min])
            if !isnan(avg_Ap_H[2])
                write(file, "<Ap>_L = $(avg_Ap_L[1]) ± $(avg_Ap_L[2])\n")
                write(file, "<Ap>_H = $(avg_Ap_H[1]) ± $(avg_Ap_H[2])\n\n")
            else
                write(file, "<Ap>_L = $(avg_Ap_L[1])\n")
                write(file, "<Ap>_H = $(avg_Ap_H[1])\n\n")
            end
        end
        if secondary_output_Tₖ == "YES"
            for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
                Tₖ_A_Z_TKE = DataFrame(
                    A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Value = Raw_output_datafile.Tₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
                )
                avg_Tₖ_L = Average_value(Tₖ_A_Z_TKE, y_A_Z_TKE, A_L_range)
                avg_Tₖ_H = Average_value(Tₖ_A_Z_TKE, y_A_Z_TKE, A_H_range)
                avg_Tₖ = Average_value(Tₖ_A_Z_TKE, y_A_Z_TKE, A_range)
                write(file, "<T_$(k)>_L = $avg_Tₖ_L\n<T_$(k)>_H = $avg_Tₖ_H\n<T_$(k)> = $avg_Tₖ\n\n")
            end
        end    
        if secondary_output_avg_εₖ == "YES"
            for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
                avg_εₖ_A_Z_TKE = DataFrame(
                    A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Value = Raw_output_datafile.Avg_εₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
                )
                avg_εₖ_L = Average_value(avg_εₖ_A_Z_TKE, y_A_Z_TKE, A_L_range)
                avg_εₖ_H = Average_value(avg_εₖ_A_Z_TKE, y_A_Z_TKE, A_H_range)
                avg_εₖ = Average_value(avg_εₖ_A_Z_TKE, y_A_Z_TKE, A_range)
                write(file, "avg_<ε_$(k)>_L = $avg_εₖ_L\navg_<ε_$(k)>_H = $avg_εₖ_H\navg_<ε_$(k)> = $avg_εₖ\n\n")
            end
        end
        if secondary_output_Eᵣ == "YES"
            for k in 1:maximum(Raw_output_datafile.No_Sequence[Raw_output_datafile.Tₖ .>= 0])
                Tₖ = Raw_output_datafile.Tₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
                aₖ = Raw_output_datafile.aₖ[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)]
                Eᵣ_A_Z_TKE = DataFrame(
                    A = Raw_output_datafile.A[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Z = Raw_output_datafile.Z[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    TKE = Raw_output_datafile.TKE[(Raw_output_datafile.No_Sequence .== k) .& (Raw_output_datafile.Tₖ .>= 0)],
                    Value = Energy_FermiGas.(aₖ, Tₖ)
                )
                avg_Eᵣ_L = Average_value(Eᵣ_A_Z_TKE, y_A_Z_TKE, A_L_range)
                avg_Eᵣ_H = Average_value(Eᵣ_A_Z_TKE, y_A_Z_TKE, A_H_range)
                avg_Eᵣ = Average_value(Eᵣ_A_Z_TKE, y_A_Z_TKE, A_range)
                write(file, "<Eᵣ_$(k)>_L = $avg_Eᵣ_L\n<Eᵣ_$(k)>_H = $avg_Eᵣ_H\n<Eᵣ_$(k)> = $avg_Eᵣ\n\n")
            end
        end
    end
    if secondary_output_T == "YES"
        avg_T_L = Average_value(T_A_Z_TKE, y_A_Z_TKE, A_L_range)
        avg_T_H = Average_value(T_A_Z_TKE, y_A_Z_TKE, A_H_range)
        avg_T = Average_value(T_A_Z_TKE, y_A_Z_TKE, A_range)
        write(file, "<T>_L = $avg_T_L\n<T>_H = $avg_T_H\n<T> = $avg_T\n\n")
    end    
    if secondary_output_avg_ε == "YES"
        avg_ε_L = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_L_range)
        avg_ε_H = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_H_range)
        avg_ε = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_range)
        write(file, "<ε>_L = $avg_ε_L\n<ε>_H = $avg_ε_H\n<ε> = $avg_ε\n\n")
    end  
    if secondary_output_TXE_Q == "YES"
        avg_Q = Average_value(Q_AH, y_A, A_H_range)
        avg_TXE = Average_value(txe_AH, y_A, A_H_range)
        write(file, "<Q> = $avg_Q\n<TXE>= $avg_TXE\n\n")
    end
    if secondary_output_E_excitation == "YES"
        avg_E_exi = Average_value(E_excitation_A, y_A, A_range)
        avg_E_exi_L = Average_value(E_excitation_A, y_A, A_L_range)
        avg_E_exi_H = Average_value(E_excitation_A, y_A, A_H_range)
        write(file, "<E*>_L = $avg_E_exi_L\n<E*>_H = $avg_E_exi_H\n<E*> = $avg_E_exi\n\n")
    end
end