#=
This part of the program is optional. It reads the output of the main program and processes it so that
it can be compared with experimental data (output data is averaged over experimental yield distributions)
=#
#####
#Compute Y(A,Z,TKE) form experimental Y(A,TKE) data
function Process_yield_data(A_0, fragmdomain::Distribution, dY::DataFrame)
    y_A_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    #Completes fragmentation domain in case data provided only for HF
    if !isassigned(dY.A[dY.A .< A_0/2], 1)
        for A in unique(sort(dY.A, rev=true))
            for Z in fragmdomain.Z[fragmdomain.A .== A_0 - A]
                P_A_Z = fragmdomain.Value[(fragmdomain.A .== A_0 - A) .& (fragmdomain.Z .== Z)][1]
                for TKE in unique(sort(dY.TKE[dY.A .== A]))
                    val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    if A != A_0 - A
                        push!(y_A_Z_TKE.A, A_0 - A)
                        push!(y_A_Z_TKE.Z, Z)
                        push!(y_A_Z_TKE.TKE, TKE)
                        push!(y_A_Z_TKE.Value, val)
                        push!(y_A_Z_TKE.σ, σ)
                    end
                end
            end
        end
    end
    for A in unique(sort(dY.A))
        for Z in fragmdomain.Z[fragmdomain.A .== A]
            P_A_Z = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
            for TKE in unique(sort(dY.TKE[dY.A .== A]))
                val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                if A != A_0 - A
                    push!(y_A_Z_TKE.A, A)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.Value, val)
                    push!(y_A_Z_TKE.σ, σ)
                else
                    push!(y_A_Z_TKE.A, A)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Value, 2*val)
                    push!(y_A_Z_TKE.σ, σ*sqrt(2))
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
function Average_over_Z(q_A_Z_TKE::Distribution, fragmdomain::Distribution)
    q_A_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in unique(q_A_Z_TKE.A)
        for TKE in unique(q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A)])
            Denominator = 0.0
            Numerator = 0.0
            for Z in q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    P_A_Z = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
                    Denominator += P_A_Z
                    Numerator += P_A_Z * value
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
function Average_over_TKE_Z(q_A_Z_TKE::Distribution, y_A_Z_TKE::Distribution)
    q_A = Distribution_unidym(Int[], Float64[], Float64[])
    for A in unique(q_A_Z_TKE.A)
        Denominator = 0.0
        Numerator = 0.0
        for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A)])
            for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE * value
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
function Average_over_A_Z(q_A_Z_TKE::Distribution, y_A_Z_TKE::Distribution)
    q_TKE = Distribution_unidym(Float64[], Float64[], Float64[])
    for TKE in unique(q_A_Z_TKE.TKE)
        Denominator = 0.0
        Numerator = 0.0
        for A in unique(q_A_Z_TKE.A[(q_A_Z_TKE.TKE .== TKE)])
            for Z in q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE * value
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
function Average_value(q_A_Z_TKE::Distribution, y_A_Z_TKE::Distribution)
    Denominator = 0.0
    Numerator = 0.0
    for A in unique(q_A_Z_TKE.A)
        for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A)])
            for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE * value
                    end
                end
            end
        end
    end
    return Numerator/Denominator
end
#Get L-H pair value from given q(A) distribution and A_H 
function Pair_value(q_A::Distribution_unidym, A_0, A_H)
    val_L = q_A.Value[q_A.Argument .== A_0-A_H][1]
    val_H = q_A.Value[q_A.Argument .== A_H][1]
    return val_L + val_H
end
#Obtain Y(Aₚ) distribution from Y(A,Z,TKE) & ν(A)
function Yield_post_neutron(y_A_Z_TKE::Distribution, ν_A::Distribution_unidym)
    y_Aₚ = Distribution_unidym(Int[], Float64[], Float64[])
    for A in ν_A.Argument
        Y_A = sum(y_A_Z_TKE.Value[y_A_Z_TKE.A .== A])
        σ_Y_A = sqrt(sum(y_A_Z_TKE.σ[y_A_Z_TKE.A .== A].^2))
        Aₚ = round(A - ν_A.Value[ν_A.Argument .== A][1])
        if !isassigned(y_Aₚ.Value[y_Aₚ.Argument .== Aₚ], 1)
            push!(y_Aₚ.Argument, Aₚ)
            push!(y_Aₚ.Value, Y_A)
            push!(y_Aₚ.σ, σ_Y_A)
        else
            y_Aₚ.Value[y_Aₚ.Argument .== Aₚ] .+= Y_A
            y_Aₚ.σ[y_Aₚ.Argument .== Aₚ] .= sqrt(sum(y_Aₚ.σ[y_Aₚ.Argument .== Aₚ].^2) + σ_Y_A^2)
        end
    end
    return y_Aₚ
end
#Obtain Y(Z,Aₚ) distribution from Y(A,Z,TKE) & ν(A,Z)
function Yield_post_neutron(y_A_Z_TKE::Distribution, ν_A_Z_TKE::Distribution)
    y_Z_Aₚ = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in unique(ν_A_Z_TKE.A)
        for Z in unique(ν_A_Z_TKE.Z[ν_A_Z_TKE.A .== A])
            Denominator = 0.0
            Numerator = 0.0
            for TKE in ν_A_Z_TKE.TKE[(ν_A_Z_TKE.A .== A) .& (ν_A_Z_TKE.Z .== Z)]
                value = ν_A_Z_TKE.Value[(ν_A_Z_TKE.A .== A) .& (ν_A_Z_TKE.Z .== Z) .& (ν_A_Z_TKE.TKE .== TKE)][1]
                if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                    Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                    Denominator += Y_A_Z_TKE
                    Numerator += Y_A_Z_TKE * value
                end
            end
            if Denominator > 0
                ν_A_Z = Numerator/Denominator
                Y_A_Z = sum(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z)])
                σ_Y_A_Z = sqrt(sum(y_A_Z_TKE.σ[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z)].^2))
                Aₚ = round(A - ν_A_Z)
                if !isassigned(y_Z_Aₚ.Value[(y_Z_Aₚ.A .== Aₚ) .& (y_Z_Aₚ.Z .== Z)], 1)
                    push!(y_Z_Aₚ.A, Aₚ)
                    push!(y_Z_Aₚ.Z, Z)
                    push!(y_Z_Aₚ.Value, Y_A_Z)
                    push!(y_Z_Aₚ.σ, σ_Y_A_Z)
                else
                    y_Z_Aₚ.Value[(y_Z_Aₚ.A .== Aₚ) .& (y_Z_Aₚ.Z .== Z)] .+= Y_A_Z
                    y_Z_Aₚ.σ[(y_Z_Aₚ.A .== Aₚ) .& (y_Z_Aₚ.Z .== Z)] .= sqrt(sum(y_Z_Aₚ.σ[(y_Z_Aₚ.A .== Aₚ) .& (y_Z_Aₚ.Z .== Z)].^2) + σ_Y_A_Z^2)
                end
            end
        end
    end
    return y_Z_Aₚ
end